function Fout = FFTn(Fin,n,fs,fmax)

Fout=zeros(1,n);
N=length(Fin);
for k=1:floor(N/n)
    Fout = Fout + abs(fft(Fin((k-1)*n+1:k*n)));
end

if nargin>2
    vf=0:fs/n:fs-fs/n;
%     figure(fig)
    plot(vf(vf<=fmax),Fout(vf<=fmax))
%     xlabel('Frequency (Hz)','FontSize',18)
%     ylabel('Amplitude (a.u.)','FontSize',18)

end

end