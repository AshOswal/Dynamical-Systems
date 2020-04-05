% simulate coherence
clear all;close all;
p.dTsg     = 6e-3;   % ms
p.dTgs     = 6e-3;   % ms
p.dTgg     = 4e-3;   % ms
p.dTcs     = 5.5e-3;
p.dTsc     = 20e-3;  %21.5e-3;
p.dTsgpi   = 6e-3;
p.dTgpegpi = 6e-3;
p.dTgpigpi = 4*1e-3;

p.Ts   = 12.8*1e-3;  % ms
p.Tg   = 20*1e-3;    % ms
p.Te   = 11e-3;      % ms - we usually use 11 here
p.Ti   = 11e-3;      % ms - we usually use 11 here
p.Tgpi = 14e-3;      % see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2585404/#!po=67.1053 for reference

p.Ctx  = 172.18;     %27;       % spm/s
p.Str  = 8.46;       %  2;      % sp/s
p.Ms   = 300;        % sp/s
p.Bs   = 10;         %17;     % sp/s
p.Mg   = 400;        % sp/s
p.Bg   = 20;         %75;     % sp/s
p.Mgi  = 300;
p.Bgi  = 18;
p.Mi   = 205;        % normal range 200-330 sp/s
p.Bi   = 9.87;       % normal range 0-20 sp/s
p.Me   = 75;         % normal range 50-80 sp/s
p.Be   = 17.85;      % normal range 0-20 sp/s

% for the cortical populations
p.dTcs  = 5.5e-3;   % ms - allowed range 1-10 ms
p.dTcc  = 6e-3;     %5e-3%1e-3;%4.65e-3;  % 2;     % ms - allowed range 1-10 ms

p.Wsg  = 4.87;     %19 + p.K*(20-19);
p.Wgs  = 1.33;     %1.12 + p.K*(10.7-1.12);
p.Wgg  = 0.53;     %0.53;%6.6 + p.K*(12.3 - 6.6);
p.Wcs  = 7.98;     %10.98;%2.42 + p.K*(9.2 - 2.42);

% cortical parameters
p.Wsc     = 0.5;   %0.4;%0.4;%0.4%8.93;
p.Wcc     = 4;
p.Wsgpi   = 2.7;%2.7
p.Wgpegpi = 1;

p.S       = 0.8;%0.8;%5.5;
p.SN      = 300;%250;% normally we use 340;%105;%405
p.PhaseP  = 0;
p.Method  = 'S';
p.Plot    = 1;
cplot     = 0;
%% show the effect of noise on coherence between cortex and STN/GPI

p.Wcs     = 8;%8;
p.Wsgpi   = 3;%2.9;%4.8;
p.length  = 40;
%p.Wgs     = 1.3;
%p.Wsg     = p.Wsg;
% try
%     load('cohs.mat');
%     load('pows.mat');
% catch
for n = 1:50
    for L = 1:3
        for ES = 1:12
            
            %p.Wsc = 0 + 0.2*(L-1);
            p.Wsc = 0 + 0.5*(L-1);
            p.Wcs = 0 + (ES-1);
            % simulate the data
            [Y3,F3,S] = simulate_PD(p);close gcf;
            cfg = [];
            timelock = [];
            
            cfg.label    = {'STN';'GPe';'GPi';'E';'I'};
            
            Tl = 2;
            nt = floor(p.length/Tl);
            
            for i = 2:nt
                sl = floor(length(Y3.x)/nt);
                cfg.time{i-1} = Y3.x(1:sl);
                cfg.trial{i-1} = Y3.y(:,sl*(i-1)+1:sl*i);
                
                %             cfg.trial{i} = randn(size(cfg.trial{i}));
            end
            timelock     = cfg;
            %figure;plot(cfg.time{1},cfg.trial{1}')
            
            cfg = [];
            cfg.output       = 'powandcsd';
            cfg.channelcmb   = {timelock.label{1} timelock.label{4};...
                timelock.label{3} timelock.label{4};...
                timelock.label{1} timelock.label{5};...
                timelock.label{3} timelock.label{5}};
            cfg.keeptrials   = 'no';
            cfg.keeptapers   = 'no';
            cfg.taper        = 'dpss';
            cfg.method       = 'mtmfft';
            cfg.pad          = 'nextpow2';
            cfg.foilim       = [2 50];
            cfg.tapsmofrq    = 5;
            inp = ft_freqanalysis(cfg, timelock);
            
            cfg.method = 'coh';
            coh = ft_connectivityanalysis(cfg,inp);
            if cplot
                figure;
                subplot(2,1,1);
                plot(coh.freq,mean(coh.cohspctrm([1 ],:),1),'b');hold on;
                plot(coh.freq,mean(coh.cohspctrm([2 ],:),1),'r');
                subplot(2,1,2);
                plot(inp.freq,inp.powspctrm([1 2 4],:));
                
            end
            cohs{ES,L,n} = coh;
            pows{ES,L,n} = inp;
            
            %{
            % granger directionality
            D = spm_eeg_ft2spm(timelock, 'temp');
            Dc = granger_direction(D, cfg.channelcmb);
            G  = Dc(:,:,:,Dc.indtrial('granger_orig')) - Dc(:,:,:,Dc.indtrial('granger_reversed'));
            
            %%
            Gall(ES,L,:,:,n)= G;
            %}
        end
    end
    disp(num2str(n));
end
save('cohs.mat','cohs');
save('pows.mat','pows');
%save('Gall.mat','Gall');

%end

    
%end
%%
for L = 1:3
    for ES = 1:12
        coh           = collect(cohs(ES,L,:),'cohspctrm');
        hb            = find(cohs{1,1,1}.freq > 22 &  cohs{1,1,1}.freq < 29);
        coh_stn_M(L,ES) = mean(mean(mean(coh(:,[1 3],hb),1),2),3);
        coh_stn_S(L,ES) = mean(mean(std(coh(:,[1 3],hb),1),2),3)./sqrt(size(coh,1));
    end
end
figure;

plot(0:11,coh_stn_M(1,:),'bo','MarkerFaceColor','b');hold on;
%plot(0:11,coh_stn_M(2,:),'bo','MarkerFaceColor','b');hold on;
%plot(0:11,coh_stn_M(3,:),'go','MarkerFaceColor','g');hold on;
% h = legend('Wgie = 0','Wgie = 0.5','Wgie = 1');
%set(h,'box','off','Location','NorthEastOutside');

eb = errorbar(0:11,coh_stn_M(1,:),coh_stn_S(1,:),'LineStyle', 'none');hold on
set(eb, 'color', 'b', 'LineWidth', 1);hold on;

%eb = errorbar(0:11,coh_stn_M(2,:),coh_stn_S(2,:),'LineStyle', 'none');hold on
%set(eb, 'color', 'b', 'LineWidth', 1);hold on;

%eb = errorbar(0:11,coh_stn_M(3,:),coh_stn_S(3,:),'LineStyle', 'none');hold on
%set(eb, 'color', 'g', 'LineWidth', 1);hold on;
set(gca,'FontSize',14);
xlim([-0.5 11.5])


%%

for L = 1:3
    % coherence
    ind_ES_7 = 9;
    coh = collect(cohs(ind_ES_7,L,:),'cohspctrm');
    coh_stn_m = squeeze(mean(mean(coh(:,[1 3],:),1),2));
    coh_stn_s = squeeze(mean(std(coh(:,[1 3],:),1),2)) ./sqrt(size(coh,1));
    
    coh_gp_m  = squeeze(mean(mean(coh(:,[2 4],:),1),2));
    coh_gp_s  = squeeze(mean(std(coh(:,[2 4],:),1),2)) ./sqrt(size(coh,1));
    
    % power
    pow   = sqrt(collect(pows(ind_ES_7,L,:),'powspctrm'));
    powSTN_m = squeeze(mean(pow(:,1,:),1));
    powSTN_s = squeeze(std(pow(:,1,:),1)./sqrt(size(pow,1)));
    
    powGPi_m = squeeze(mean(pow(:,3,:),1));
    powGPi_s = squeeze(std(pow(:,3,:),1)./sqrt(size(pow,1)));
    
    powCTX_m = squeeze(mean(mean(pow(:,[4 5],:),1),2));%squeeze(mean(mean(squeeze(pow(:,[4 5],:)),2),1))';
    powCTX_s = squeeze(mean(std(pow(:,[4 5],:),1),2)./sqrt(size(pow,1)));
    %
    figure;
    subplot(1,2,1)
    plot(pows{1}.freq,powSTN_m,'b');hold on;
    plot(pows{1}.freq,powGPi_m,'r');hold on;
    plot(pows{1}.freq,powCTX_m,'k');hold on

    
    fill([pows{1}.freq';flipud(pows{1}.freq')],[(powSTN_m-powSTN_s);flipud((powSTN_m +powSTN_s))],'b','linestyle','none','FaceAlpha',0.3);
    fill([pows{1}.freq';flipud(pows{1}.freq')],[(powGPi_m-powGPi_s);flipud((powGPi_m +powGPi_s))],'r','linestyle','none','FaceAlpha',0.3);
    fill([pows{1}.freq';flipud(pows{1}.freq')],[(powCTX_m-powCTX_s);flipud((powCTX_m +powCTX_s))],'k','linestyle','none','FaceAlpha',0.3);
    set(gca,'FontSize',14);
    h = legend({'STN','GPi','CTX'});
    set(h,'box','off');
   
    axis([5 45 10 15])
    %xlim([5 45])
    xlabel('Frequency (Hz)');
    ylabel('Power (AU)');
    set(gca,'FontSize',14,'XTick',0:10:45);
    subplot(1,2,2)
    f = cohs{2}.freq;
    s1= coh_stn_s;
    plot(f,coh_stn_m,'b');hold on;
    plot(f,coh_gp_m,'r');hold on
   
    fill([f';flipud(f')],[coh_stn_m-coh_stn_s;flipud(coh_stn_m+coh_stn_s)],'b','linestyle','none','FaceAlpha',0.3);
    fill([f';flipud(f')],[coh_gp_m-coh_gp_s;flipud(coh_gp_m+coh_gp_s)],'r','linestyle','none','FaceAlpha',0.3);
    axis([5 45 0 0.2]);
    
    xlabel('Frequency (Hz)');
    ylabel('Coherence (AU)');
    set(gca,'FontSize',14,'XTick',0:10:45);
    xlim([5 45])
    h = legend({'CTX-STN','CTX-GPi'});
    set(h,'box','off');
   
end


%%
%Gall = mean(Gall,4);
lb = find(Dc.frequencies<15 & Dc.frequencies >10);   
hb = find(Dc.frequencies>20 & Dc.frequencies <30); 
allb = find(Dc.frequencies>10 & Dc.frequencies <30);
%%
for L = 1:6
    % mean granger spectra
    
    
    CTX_STN_HBm = mean(mean(mean(Gall(L,[2,6],hb,:),2),3),4);
    CTX_STN_HBs = std(squeeze(mean(mean(Gall(L,[2,6],hb,:),2),3)),1)/sqrt(20);
    
    CTX_STN_LBm = mean(mean(mean(Gall(L,[2,6],lb,:),2),3),4);
    CTX_STN_LBs = std(squeeze(mean(mean(Gall(L,[2,6],lb,:),2),3)),1)/sqrt(20);
    
    STN_CTX_HBm = mean(mean(mean(Gall(L,[1,5],hb,:),2),3),4);
    STN_CTX_HBs = std(squeeze(mean(mean(Gall(L,[1],hb,:),2),3)),1)/sqrt(20);
   
    STN_CTX_LBm = mean(mean(mean(Gall(L,[1,5],lb,:),2),3),4);
    STN_CTX_LBs = std(squeeze(mean(mean(Gall(L,[1],lb,:),2),3)),1)/sqrt(20);
  
    CTX_GPi_m = mean(mean(mean(Gall(L,[4],allb,:),2),3),4);
    CTX_GPi_s = std(squeeze(mean(mean(Gall(L,[4],allb,:),2),3)),1)/sqrt(20);

    GPi_CTX_m = mean(mean(mean(Gall(L,[3],allb,:),2),3),4);
    GPi_CTX_s = std(squeeze(mean(mean(Gall(L,[3],allb,:),2),3)),1)/sqrt(20);
    
    means1 = [CTX_STN_HBm,STN_CTX_HBm;CTX_STN_LBm,STN_CTX_LBm];
    means2 = [CTX_GPi_m,GPi_CTX_m;0,0];
    st     = [CTX_STN_HBs,STN_CTX_HBs;CTX_STN_LBs,STN_CTX_LBs];
    st2    = [CTX_GPi_s,GPi_CTX_s;0,0];
    
    figure;
    hBar = barwitherr(st,[1 2 ],means1,'EdgeColor',[0 0 0]);
    set(hBar(1),'facecolor',[0 0 1]);
    set(hBar(2),'facecolor',[0.5 0.5 0.7]);
    h = legend('CTX -> STN','STN -> CTX');
    set(h,'box','off','Location','NorthEastOutside');
    set(gca, 'FontSize',12,'XTickLabel',[]);

    figure;
    hBar3 = barwitherr(st2,[1 2 ],means2,'EdgeColor',[0 0 0]);
    set(hBar3(1),'facecolor',[1 0 0]);
    set(hBar3(2),'facecolor',[0.7 0.5 0.5]);
    h = legend('CTX -> GPi','GPi -> CTX');
    set(h,'box','off','Location','NorthEastOutside');
    set(gca, 'FontSize',12,'XTickLabel',[]);

    

 

    

end
