function Dc = granger_direction(D, channelcmb)

    
odata = D.fttimelock(sort(D.indchannel(unique(channelcmb))), ':', ':');
rdata = odata;
rdata.trial = rdata.trial(:, :, end:-1:1);
sdata = odata;

ind = spm_match_str(sdata.label, sdata.label);

sdata.trial(:, ind, :) = sdata.trial([2:end 1], ind, :);
data = {odata, rdata, sdata};
%%
for i = 1:numel(data)
    
    fstep = 1/(D.nsamples/D.fsample);
    fstep = round(1/fstep)*fstep;
    
    %fstep = 1;
    
    foi     = 0:fstep:D.fsample/2;
    foi     = foi(1:(end-1));
    fres    = 0*foi+2.5;
    fres(fres>25) = 0.1*fres(fres>25);
    fres(fres>50) = 5;
    
    cfg = [];
    cfg.output ='fourier';
    cfg.channelcmb=channelcmb;
    
    cfg.keeptrials = 'yes';
    cfg.keeptapers='yes';
    cfg.taper = 'dpss';
    cfg.method          = 'mtmfft';
    cfg.foi     = foi;
    cfg.tapsmofrq = fres;
    %cfg.pad = 20;
    
    inp{i} = ft_freqanalysis(cfg, data{i});
    %
    cfg = [];
    cfg.channelcmb=channelcmb;
    cfg.method  = 'coh';
    
    res{1, i} = ft_connectivityanalysis(cfg, inp{i});
    
    %
    cfg.complex = 'imag';
    res{2, i} = ft_connectivityanalysis(cfg, inp{i});
    
    %
    cfg.complex = 'real';
    res{3, i} = ft_connectivityanalysis(cfg, inp{i});
    
    cfg = rmfield(cfg, 'complex');
    %
    cfg.method = 'granger';
    %cfg.granger.init = 'chol';
    cfg.granger.sfmethod ='bivariate';  
    
    res{4, i} = ft_connectivityanalysis(cfg, inp{i});
    
    cfg.method = 'instantaneous_causality';
    res{5, i} = ft_connectivityanalysis(cfg, inp{i});
end
%%
Nchannels = 2*size(channelcmb, 1);
Nfrequencies = length(res{1, 1}.freq);
Ntrials = numel(res);

Dc = clone(D, ['C' D.fname], [Nchannels Nfrequencies 1 Ntrials]);
Dc = Dc.frequencies(':', res{1, 1}.freq);
Dc = timeonset(Dc, 0);
Dc = fsample(Dc, 1);
Dc = transformtype(Dc, 'TF');

reslabels = {'coh', 'coh', 'coh', 'granger', 'instant'};
outlabels = {'coh', 'imagcoh', 'realcoh', 'granger', 'instant'};
datalabels = {'orig', 'reversed', 'shifted'};
cl = {};
chanl = {};

for i = 1:numel(reslabels)
    for j = 1:numel(datalabels)
        trialind = sub2ind([numel(reslabels), numel(datalabels)], i, j);
        cl{trialind, 1} = [outlabels{i} '_' datalabels{j}];
        for k = 1:size(channelcmb, 1)            
            ind1 = intersect(strmatch(channelcmb{k, 1}, res{i, j}.labelcmb(:, 1)), ...
                strmatch(channelcmb{k, 2}, res{i, j}.labelcmb(:, 2)));
            ind2 = intersect(strmatch(channelcmb{k, 2}, res{i, j}.labelcmb(:, 1)), ...
                strmatch(channelcmb{k, 1}, res{i, j}.labelcmb(:, 2)));
            if isequal(reslabels{i}, 'coh') 
                ind1 = [ind1 ind2];
                ind2 = ind1;
            end                        
            
            Dc(2*k-1, :, :, trialind) = res{i, j}.([reslabels{i} 'spctrm'])(ind1, :);
            Dc(2*k, :, :, trialind)   = res{i, j}.([reslabels{i} 'spctrm'])(ind2, :);
            
            chanl{2*k-1, 1} = [channelcmb{k, 1} '->' channelcmb{k, 2}];
            chanl{2*k,   1} = [channelcmb{k, 2} '->' channelcmb{k, 1}];
        end
        
    end
end
Dc = chanlabels(Dc, ':', chanl);
Dc = conditions(Dc, ':', cl);

save(Dc);
%%

return
%%

ROI = {'CTX_theta_LFP_R23'};


cnd = {'granger_orig', 'granger_reversed',  'granger_shifted'};%  Dc.condlist;%{'coh_orig', 'coh_shifted'};
%cnd = {'instant_orig', 'instant_reversed',  'instant_shifted'};
%cnd = {'coh_orig', 'coh_shifted'};
%cnd = {'imagcoh_orig', 'imagcoh_shifted'}
%cnd = {'realcoh_orig', 'realcoh_shifted'}

spm_figure('GetWin', [initials druglbl{drug+1} '_' cnd{1}]);clf;

trialind = Dc.indtrial(cnd);

clf;
for i = 1:size(ROI, 1)
    ind1 = strmatch(ROI{i, 1}, Dc.chanlabels);
    ind2 = strmatch(fliplr(ROI{i, 1}), cellfun(@fliplr, Dc.chanlabels',  'UniformOutput', false));
    
    subplot(size(ROI, 1), 2, 2*i-1);
    plot(Dc.frequencies, squeeze(mean(Dc(ind1, :, :, trialind), 1)));
    
    xlim([5 45]);
    
    if i == 1
        legend(Dc.conditions(trialind), 'Interpreter', 'none')
    end
    
    title([ROI{i, 1} '->LFP']);
    
    subplot(size(ROI, 1), 2, 2*i);
    plot(Dc.frequencies, squeeze(mean(Dc(ind2, :, :, trialind), 1)));  
    
    xlim([5 45]);
    title(['LFP->' ROI{i, 1}]);    
end
%%
keyboard
