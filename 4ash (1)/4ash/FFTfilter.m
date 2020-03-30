function DataF=FFTfilter(Data,Fs,Fhp,Flp)

N = length(Data);

dF  = Fs/N;
vF  = 0:dF:floor(Fs/2);
vF  = [vF -vF(ceil(N/2):-1:2)];
ind = logical((abs(vF)<Fhp)+(abs(vF)>Flp));

FData = fft(Data,[],2);
FData(:,ind) = 0;
DataF = ifft(FData,[],2);

end