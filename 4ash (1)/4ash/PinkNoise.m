function PN = PinkNoise(N,b)

Wn = randn(N,1);
if nargin>1
    b = .5*b;
else
    b = .5;
end
Pc = (1./((.25:N).^(b)))';
Pc(N:-1:(N-round(N/2))+2) = Pc(2:round(N/2));
Pc = Pc/norm(Pc);
PN = real(ifft(fft(Wn).*(Pc)));

end