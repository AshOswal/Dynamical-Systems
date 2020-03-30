function []=FFTplot(Vin,fs,a,b,str,c)
% FFTplot(Vin,fs,~)

N=length(Vin);
vf = 0:fs/N:fs-fs/N;

if nargin ==6
    [VinF vf]= pmtm(Vin,4);
    VinF = sqrt(VinF);
    vf = vf*fs/2/pi;
else
    VinF=fft(Vin);
end

if nargin<5
    str='b';
end
if nargin>3 && b
    plot(log(vf(1:ceil(N/2))),log(abs(VinF(1:ceil(N/2)))),str)
    xlabel('Log Frequency (Hz)','FontSize',14)
    ylabel('Log Amplitude (a.u.)','FontSize',14)
elseif nargin>2 && a
    plot(vf(1:ceil(N/2)),log(abs(VinF(1:ceil(N/2)))),str)
    xlabel('Frequency (Hz)','FontSize',14)
    ylabel('Log Amplitude (a.u.)','FontSize',14)
else
    plot(vf(1:ceil(N/2)),abs(VinF(1:ceil(N/2))),str)
    xlabel('Frequency (Hz)','FontSize',14)
    ylabel('Amplitude (a.u.)','FontSize',14)
end

end