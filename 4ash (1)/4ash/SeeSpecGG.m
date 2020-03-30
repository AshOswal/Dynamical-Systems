
impF = [ 0.5 4 4 8 8 12 12 30 30 50 50 100 ]'; % Frequency limmit pairs for the frequency bands of intrest

% Parameters for TF analysis

R = 7; % Morlet Wavelet factor
N = 50; % # Frequency Samples
maxf = impF(end); % S.fsample / 2 or other;
minf = impF(1)-0.1; % Defined by me
alfa = (maxf/minf)^(1/N)-1; % According to the expression achived by fn = ((1+1/R)^n)*f0 where 1/R = alfa

linear_f_step = 2;

vf = ((1+alfa).^(1:N)) * minf;
S.fsample=250;


%%

% for iC=1:64;
for iC=[1:9];
    load([num2str(iC) '.mat'])
    Dtf_f=log(abs(S3)+0.0001);
    % close all
     figure('Name',['Componente ' num2str(iC)],'Position',[10, 550, 1900, 400],...
        'NumberTitle','off')
    set(gca,'Units', 'pixels','Position', [20,20,1870,370])
    imagesc(squeeze(Dtf_f(:,:)));
    Fplot = (log([1 2 4 8 16 32 64])-log(vf(1)))./log(1+alfa)+1;
    N = length(Dtf_f);
    aux0 = N:-1:1;
    aux1=1:N;
    
    hold on
    plot([1 N],[Fplot',Fplot'],'k')
    for k=1:length(GIlist)
        plot([EEG.event(GIlist(k)).latency EEG.event(GIlist(k)).latency;],[0 100],'r','LineWidth',2)
    end
    for k=1:length(EIlist)
        plot([EEG.event(EIlist(k)).latency EEG.event(EIlist(k)).latency;],[0 100],'k','LineWidth',2)
    end
    for k=1:length(LIlist)
        plot([EEG.event(LIlist(k)).latency EEG.event(LIlist(k)).latency;],[0 100],'w','LineWidth',2)
    end
    for k=1:length(HeadLsit)
        plot([EEG.event(HeadLsit(k)).latency EEG.event(HeadLsit(k)).latency;],[0 100],'b','LineWidth',2)
    end
    for k=1:length(SlowList)
        plot([EEG.event(SlowList(k)).latency EEG.event(SlowList(k)).latency;],[0 100],'r','LineWidth',2)
    end
    for k=1:length(SlowList)
        plot([EEG.event(SlowList(k)).latency+EEG.event(SlowList(k)).duration, EEG.event(SlowList(k)).latency+EEG.event(SlowList(k)).duration;],[0 100],'w','LineWidth',2)
    end
    for k=1:length(ScanList)
        plot([EEG.event(ScanList(k)).latency EEG.event(ScanList(k)).latency;],[0 100],'k','LineWidth',1)
    end
%     figure
    set(gca,'YTick',Fplot)
    set(gca,'YTickLabel',[1 2 4 8 16 32 64])
    set(gca,'XTick',(0:5000:N))
    set(gca,'XTickLabel',(0:5000:N)/S.fsample)
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    plot((log(max(FFTfilter(RegRMSF,S.fsample,0,0.5),vf(1)))-log(vf(1)))./log(1+alfa)+1,'k','LineWidth',2)
    hold off
end


%%

iC = 9;
load([num2str(iC) '.mat'])

GIlist=[];
EIlist=[];
LIlist=[];
HeadLsit=[];
ScanList=[];
ScanTimes=[];
SlowList=[];
RSList=[];
RSIList=[];


for k=1:length(EEG.event);
    if strcmp(EEG.event(k).type,'GI') GIlist=[GIlist; k]; end
    if strcmp(EEG.event(k).type,'gelastic') GIlist=[GIlist; k]; end
    if strcmp(EEG.event(k).type,'EI_LF') EIlist=[EIlist; k]; end
    if strcmp(EEG.event(k).type,'LI_LF') LIlist=[LIlist; k]; end
    if strcmp(EEG.event(k).type,'EI') EIlist=[EIlist; k]; end
    if strcmp(EEG.event(k).type,'LI') LIlist=[LIlist; k]; end
    if strcmp(EEG.event(k).type,'rhythmic_slow') RSList=[RSList; k]; end
    if strcmp(EEG.event(k).type,'Interictal') RSList=[RSList; k]; end
    if strcmp(EEG.event(k).type,'rhythmic_slow_intermixed_1-2HZ') RSIList=[RSIList; k]; end
    if strcmp(EEG.event(k).type,'head') HeadLsit=[HeadLsit; k]; end
    if strcmp(EEG.event(k).type,'slow') SlowList=[SlowList; k]; end
    
    if strcmp(EEG.event(k).type,'Scan Start')
        ScanList=[ScanList; k];
        ScanTimes=[ScanTimes, EEG.event(k).latency ];
    end
end

ScanTimes=ScanTimes+mean(diff(ScanTimes))/2;
ScanTimes = round([ScanTimes(1)-mean(diff(ScanTimes)), ScanTimes]);
% ScanTimes = round([ScanTimes]);

GIReg=zeros(1,length(S3));
EIReg=zeros(1,length(S3));
LIReg=zeros(1,length(S3));
RSReg=zeros(1,length(S3));
RSIReg=zeros(1,length(S3));
for k=1:length(GIlist)
    GIReg(EEG.event(GIlist(k)).latency:EEG.event(GIlist(k)).latency+EEG.event(GIlist(k)).duration)=1;
end
for k=1:length(EIlist)
    EIReg(EEG.event(EIlist(k)).latency:EEG.event(EIlist(k)).latency+EEG.event(EIlist(k)).duration)=1;
end
for k=1:length(LIlist)
    LIReg(EEG.event(LIlist(k)).latency:EEG.event(LIlist(k)).latency+EEG.event(LIlist(k)).duration)=1;
end
for k=1:length(RSList)
    RSReg(EEG.event(RSList(k)).latency:EEG.event(RSList(k)).latency+EEG.event(RSList(k)).duration)=1;
end
for k=1:length(RSIList)
    RSIReg(EEG.event(RSIList(k)).latency:EEG.event(RSIList(k)).latency+EEG.event(RSIList(k)).duration)=1;
end

[ hrf , p ] = spm_hrf (1/S.fsample); hrftd = diff(hrf);

HGIReg = conv(hrf,GIReg);
HEIReg = conv(hrf,EIReg);
HLIReg = conv(hrf,LIReg);
HRSReg = conv(hrf,RSReg);
HRSIReg = conv(hrf,RSIReg);
HGIReg = HGIReg-mean(HGIReg); HGIReg=HGIReg/norm(HGIReg);
HGIReg_DS = HGIReg(ScanTimes);
HLIReg = HLIReg-mean(HLIReg); HLIReg=HLIReg/norm(HLIReg);
HLIReg_DS = HLIReg(ScanTimes);
HEIReg = HEIReg-mean(HEIReg); HEIReg=HEIReg/norm(HEIReg);
HEIReg_DS = HEIReg(ScanTimes);
HRSReg = HRSReg-mean(HRSReg); HRSReg=HRSReg/norm(HRSReg);
HRSReg_DS = HRSReg(ScanTimes);
HRSIReg = HRSIReg-mean(HRSIReg); HRSIReg=HRSIReg/norm(HRSIReg);
HRSIReg_DS = HRSIReg(ScanTimes);

figure;plot([HEIReg_DS' HGIReg_DS' HLIReg_DS' HRSReg_DS' HRSIReg_DS'])

fid = fopen(['GIRegBF1.txt'],'w');
fprintf(fid,'%d\n' ,HGIReg_DS');
fclose(fid);
fid = fopen(['EIRegBF1.txt'],'w');
fprintf(fid,'%d\n' ,HEIReg_DS');
fclose(fid);
fid = fopen(['LIRegBF1.txt'],'w');
fprintf(fid,'%d\n' ,HLIReg_DS');
fclose(fid);
fid = fopen(['RSRegBF1.txt'],'w');
fprintf(fid,'%d\n' ,HRSReg_DS');
fclose(fid);
fid = fopen(['RSIRegBF1.txt'],'w');
fprintf(fid,'%d\n' ,HRSIReg_DS');
fclose(fid);


%% Extract Metrics

iC = 9;

load([num2str(iC) '.mat'])
% Frequ?ncia * Amplitude -------------
% RegFP = zeros(length(compintrest),length(S3));
% Pot?ncia Total --------------------------
RegTP = zeros(1,length(S3));
% Frequ?ncia RMS -------------
RegRMSF = zeros(1,length(S3));
RegRMSFTP= zeros(1,length(S3));
% Frequ?ncia uRMS -------------
% ReguRMSF = zeros(length(compintrest),length(S3));
% Frequencia m?dia ------------------------
% RegMF = zeros(length(compintrest),length(S3));

Fvec = vf;
FL = diff(Fvec);
FLvec = [FL/2 0] + [0 FL/2];

% RegFP(n,:) = abs(S3)'*(Fvec'.*FLvec');
RegTP(:) = abs(S3).^2'*FLvec';
RegRMSFTP(:) = sqrt( abs(S3).^2'*(Fvec.^2'.*FLvec'));
RegRMSF(:) = sqrt( abs(S3).^2'*(Fvec.^2'.*FLvec')./RegTP(:) );

% ReguRMSF(n,:) = sqrt( abs(S3).^2'*(Fvec.^2'.*FLvec') );
% RegMF(n,:) = abs(S3).^2'*(Fvec'.*FLvec')./RegTP(n,:)';

[ hrf , p ] = spm_hrf (1/S.fsample); hrftd = diff(hrf);
hrf=hrf/norm(hrf);

HRMSF   = conv(hrf,RegRMSF);    % Principal
HRMSF=HRMSF-mean(HRMSF); HRMSF=HRMSF/norm(HRMSF);
HRMSFTP   = conv(hrf,RegRMSFTP);    % Principal
HRMSFTP=HRMSFTP-mean(HRMSFTP); HRMSFTP=HRMSFTP/norm(HRMSFTP);
HTP =  conv(hrf,RegTP);
HTP=HTP-mean(HTP); HTP=HTP/norm(HTP);

HRMSF_DS= HRMSF(ScanTimes);
HRMSFTP_DS= HRMSFTP(ScanTimes);
HTP_DS = HTP(ScanTimes);

% HTDRMSF  = conv([hrftd; 0],RegRMSF);  % Temporal Derivative
% HTDRMSF=HTDRMSF-mean(HTDRMSF); HTDRMSF=HTDRMSF/norm(HTDRMSF);
% HTDRMSF_DS= HTDRMSF(ScanTimes);


figure;plot([HRMSF_DS' HRMSFTP_DS' HTP_DS'])

fid = fopen([num2str(iC) 'RegBF1.txt'],'w');
fprintf(fid,'%d\n' ,HRMSF_DS');
fclose(fid);
fid = fopen(['RegRMSFTP' num2str(iC) 'RegBF1.txt'],'w');
fprintf(fid,'%d\n' ,HRMSFTP_DS');
fclose(fid);
fid = fopen(['RegTP' num2str(iC) 'RegBF1.txt'],'w');
fprintf(fid,'%d\n' ,HTP_DS');
fclose(fid);
% fid = fopen([num2str(iC) 'RegBF2.txt'],'w');
% fprintf(fid,'%d\n' ,HTDRMSF_DS');
% fclose(fid);

%%

