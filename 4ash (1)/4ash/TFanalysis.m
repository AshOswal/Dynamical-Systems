function [S3]=TFanalysis(AD,R,Nf,minf,maxf,fsample,arbn)
% Time Frequency analysis
% TFanalysis(AD,R,Nf,minf,maxf,fsample,arbn)
% AD      - data in row format
% R       - Wavelet Factor                          (e.g. 7)
% Nf      - Number of frequency bins                (e.g 124)
% minf    - Minimum frequency of analysis           (e.g 1 Hz)
% maxf    - Maximum Frequency of analysis           (e.g 124 Hz)
% fsample - Data sampling Frequency                 (e.g 256 Hz)
% arbn    - Ploting soft-threshold arbitrary number (e.g. 0.001 for typical icEEG SNR values)

alfa = (maxf/minf)^(1/Nf)-1; % According to the expression achived by fn = ((1+1/R)^n)*f0 where 1/R = alfa
vf = ((1+alfa).^(1:Nf)) * minf;

dataeegF = fft(AD(1,:));
N    = length(dataeegF);
S3   = zeros(length(vf),N);
n=1;
for F=vf
    Fwavelet = (MorletWavelet(F,R,fsample,N));
    S3(n,:)  = ifft(dataeegF.*fft(Fwavelet))./N;
    n = n+1;
end

if nargin>6
    Dtf_f=log(abs(S3)+arbn);
    y=[1 Nf];
    x =[1/fsample length(AD)/fsample];
    figure('Name','TF analysis','NumberTitle','off')
    % set(gca,'Position', [0.04,.12,.95,.87]); % 'Units', 'pixels',
    imagesc(x,y,squeeze(Dtf_f(:,:)));
    Fplot = (log([1 2 4 8 16 32 64 128 256])-log(vf(1)))./log(1+alfa)+1;
    hold on
    plot(x,[Fplot',Fplot'],'k')
    hold off
    set(gca,'YTick',Fplot)
    set(gca,'YTickLabel',[1 2 4 8 16 32 64 128 256],'FontSize',16)
    % set(gca,'XTick',(0:5000:N))
    % set(gca,'XTickLabel',round((0:5000:N)/fsample),'FontSize',16)
    xlabel('Time (s)','FontSize',18)
    ylabel('Frequency (Hz)','FontSize',18)
end

end
