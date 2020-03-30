%% PrepareDataset

clear
loadlist = {'DATA/CA1_PC/Matfile';'DATA/CA1_PC/non_phase_coupled/Matfile';
            'DATA/CA1_OA_IN/Matfile';'DATA/CA1_RAD_IN/Matfile';
            'DATA/CA1_RAD_IN/non_phase_coupled/Matfile'};

mkdir('PreProcData')
savelist = {'PreProcData/PCM_CA1_PC';'PreProcData/PCM_CA1_PC_nopc';
            'PreProcData/CA1_OA_IN';'PreProcData/CA1_RAD_IN';
            'PreProcData/CA1_RAD_IN_nopc'};

for name = 3%1:length(loadlist)
    clear PC PCM
    load(loadlist{name})
    fprintf(['\n Processing ' loadlist{name} ' - (%i out of %i) \n'], name,length(loadlist))
    
for n =1:length(PC);
    
    ibhib = 0;
    excit = 0;
    fprintf('Processing neuron number %i\n', n)
    
    lfp = PC(n).I(:,3);
    exi = PC(n).I(:,2);
    tt  = PC(n).I(:,1);
    
    fs = 6000;
    % fs = 1/mean(diff(PC(3).E(:,1)));
    
    % lfp = lfp/max(lfp); exi = exi-mean(exi); exi=exi/max(exi); exi = exi+1;
    % plot([lfp, exi])
    minf = 16;maxf = 64;Nf = 200;
    [S3]=TFanalysis(lfp',4,Nf,minf,maxf,fs);
    binmax=97;
    [~, binmax] = max(sum(abs(S3),2));
    
    alfa = (maxf/minf)^(1/Nf)-1; % According to the expression achived by fn = ((1+1/R)^n)*f0 where 1/R = alfa
    vf   = ((1+alfa).^(1:Nf)) * minf;
    f0   = vf(binmax);
    
    
    lfpG            = GaussFilter(lfp',fs,4,500,10);
    InstPhase       = (angle(S3(binmax,:)));
    InstPhase       = diff(InstPhase)+pi;
    InstFrequency   = .5*([mod(InstPhase,2*pi) pi] +[pi mod(InstPhase,2*pi)])-pi;
    x               = cumsum(InstFrequency)+InstPhase(1)+pi;
    yi              = interp1(x,lfp,x(1):pi/100:x(end));
    yiE             = interp1(x,exi,x(1):pi/100:x(end));
    Lperiods = 200;
    Nperiods = floor(x(end)/2/pi)-3;
    I = Lperiods:Lperiods:Nperiods*Lperiods;
    
    Yi       = zeros(Nperiods,Lperiods);
    YiE      = zeros(Nperiods,Lperiods);
    for k = 1:Nperiods
        Yi(k,:)  = yi(I(k):I(k)+Lperiods-1);
        YiE(k,:) = yiE(I(k):I(k)+Lperiods-1);
    end
    
   
    PCM(n).IM  = YiE;
    PCM(n).ILM = Yi;
    PCM(n).If  =  f0;
    
    %----------------------------------------------------------------------
    
    lfp = PC(n).E(:,3);
    exi = PC(n).E(:,2);
    tt  = PC(n).E(:,1);
    
    fs = 6000;
    % fs = 1/mean(diff(PC(3).E(:,1)));
    
    % lfp = lfp/max(lfp); exi = exi-mean(exi); exi=exi/max(exi); exi = exi+1;
    % plot([lfp, exi])
    minf = 16;maxf = 64;Nf = 200;
    [S3]=TFanalysis(lfp',4,Nf,minf,maxf,fs);
    binmax=97;
    [~, binmax] = max(sum(abs(S3),2));
    
    alfa = (maxf/minf)^(1/Nf)-1; % According to the expression achived by fn = ((1+1/R)^n)*f0 where 1/R = alfa
    vf   = ((1+alfa).^(1:Nf)) * minf;
    f0   = vf(binmax);
    
    
    lfpG            = GaussFilter(lfp',fs,4,500,10);
    InstPhase       = (angle(S3(binmax,:)));
    InstPhase       = diff(InstPhase)+pi;
    InstFrequency   = .5*([mod(InstPhase,2*pi) pi] +[pi mod(InstPhase,2*pi)])-pi;
    x               = cumsum(InstFrequency)+InstPhase(1)+pi;
    yi              = interp1(x,lfp,x(1):pi/100:x(end));
    yiE             = interp1(x,exi,x(1):pi/100:x(end));
    Lperiods = 200;
    Nperiods = floor(x(end)/2/pi)-3;
    I = Lperiods:Lperiods:Nperiods*Lperiods;
    
    Yi       = zeros(Nperiods,Lperiods);
    YiE      = zeros(Nperiods,Lperiods);
    for k = 1:Nperiods
        Yi(k,:)  = yi(I(k):I(k)+Lperiods-1);
        YiE(k,:) = yiE(I(k):I(k)+Lperiods-1);
    end
    
   if mean(YiE(:))<mean(PCM(n).IM(:))
        PCM(n).EM  = YiE;
        PCM(n).ELM = Yi;
        PCM(n).Ef  =  f0;
   else
       fprintf('Change in E-I recordings...\n')
        PCM(n).EM  = PCM(n).IM;
        PCM(n).ELM = PCM(n).ILM;
        PCM(n).Ef  = PCM(n).If;
        PCM(n).EM  = YiE;
        PCM(n).ELM = Yi;
        PCM(n).Ef  =  f0;
   end
        
    %---------------------------------------------------------------------   
    
    lfp = PC(n).spikes(:,3);
    exi = PC(n).spikes(:,2);
    tt  = PC(n).spikes(:,1);
    
    fs = 6000;
    % fs = 1/mean(diff(PC(3).E(:,1)));
    
    % lfp = lfp/max(lfp); exi = exi-mean(exi); exi=exi/max(exi); exi = exi+1;
    % plot([lfp, exi])
    minf = 16;maxf = 64;Nf = 200;
    [S3]=TFanalysis(lfp',4,Nf,minf,maxf,fs);
    binmax=97;
    [~, binmax] = max(sum(abs(S3),2));
    
    alfa = (maxf/minf)^(1/Nf)-1; % According to the expression achived by fn = ((1+1/R)^n)*f0 where 1/R = alfa
    vf   = ((1+alfa).^(1:Nf)) * minf;
    f0   = vf(binmax);
    
    lfpG            = GaussFilter(lfp',fs,4,500,10);
    InstPhase       = (angle(S3(binmax,:)));
    InstPhase       = diff(InstPhase)+pi;
    InstFrequency   = .5*([mod(InstPhase,2*pi) pi] +[pi mod(InstPhase,2*pi)])-pi;
    x               = cumsum(InstFrequency)+InstPhase(1)+pi;
    yi              = interp1(x,lfp,x(1):pi/100:x(end));
    yiE             = interp1(x,exi,x(1):pi/100:x(end));
    Lperiods = 200;
    Nperiods = floor(x(end)/2/pi)-3;
    I = Lperiods:Lperiods:Nperiods*Lperiods;
    
    Yi       = zeros(Nperiods,Lperiods);
    YiE      = zeros(Nperiods,Lperiods);
    for k = 1:Nperiods
        Yi(k,:)  = yi(I(k):I(k)+Lperiods-1);
        YiE(k,:) = yiE(I(k):I(k)+Lperiods-1);
    end
    
    [pks I2]= findpeaks(yiE);
    
    nk = 3;
    while 1
        try
            [IDX,C,sumd,D] = kmeans(pks,nk);
            
            [C ind] = max(C);
            
            indspikes = find(IDX==ind);
            indD =  [Lperiods; diff(indspikes)];
            if sum(indD<Lperiods/10)<length(indD)/5;
                break;
            end
        end
        nk = nk+1;
        fprintf('%i K-clusters failed - trying %i K-clusters...\n',nk-1,nk);
        if nk>10; disp('K-means failed in finding spikes'); break;end
    end
%     indspikes(indD<Lperiods/10) = [];
    
    figure
    plot(yiE); hold on; plot(I2(indspikes),yiE(I2(indspikes)),'ro')
    title([loadlist{name} ' ' num2str(n)])
    drawnow
    
    yiEi = zeros(size(yiE));
    yiEi(I2(indspikes)) = 1;
    YiEi = zeros(size(YiE));
    for k = 1:Nperiods
        YiEi(k,:) = yiEi(I(k):I(k)+Lperiods-1);
    end
    
    PCM(n).SM = YiEi;
    PCM(n).SLM  = Yi;
    PCM(n).Sf  =  f0;


end


save(savelist{name},'PCM')
fprintf(['\n DONE' savelist{name} '\n'])
end





