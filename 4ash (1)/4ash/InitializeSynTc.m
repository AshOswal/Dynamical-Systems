

t = (1:length(Y.y))'/(length(Y.y)*Y.y(1,5));

T = exp(P.P(1).T); %([10,16,10,16]/1000); 
T = [10 10 10 10]/1000;

Ikernel = exp(-t/T(2))/(1-exp(-t(end)/T(2)));
Ekernel = exp(-t/T(1))/(1-exp(-t(end)/T(1)));

GI = ifft(fft(Y.y(:,7)).*fft(Ikernel))/length(Y.y);
GE = ifft(fft(Y.y(:,6)).*fft(Ekernel))/length(Y.y);

CI =  GI*(M.opt.popCurrents(1,2) - P.P(1).Vg(2));
CE =  GE*(M.opt.popCurrents(2,2) - P.P(1).Vg(1));

CI = CI - mean(CI);
CE = CE - mean(CE);

normCI = sqrt(sum(CI.^2));
normCE = sqrt(sum(CE.^2));
normY  = sqrt(sum(Y.y.^2));

M.P = M.pE;

M.P.A(1).M(1,2) = log(normY(4)/normCE);
M.P.A(2).M(2,1) = log(normY(1)/normCI);
M.P.A(3).M(2,2) = log(normY(2)/normCE);
M.P.A(4).M(1,1) = log(normY(3)/normCI);

%P.P(2).Vg  = [0,-70,0,-70];             % [mV] Reversal portential of the families of synaptic channels 

M.P.C(1).M = [log(.1) log(8)]';
M.P.C(2).M = [log(.1) log(.1)]';
M.P.C(3).M = [-32 -32]';
M.P.C(4).M = [-32 -32]';

M.opt.Drawplots=0;



YCI1 = Y.y(:,1)/sqrt(sum(Y.y(:,1).^2));
YCE1 = Y.y(:,2)/sqrt(sum(Y.y(:,2).^2));

YCI2 = Y.y(:,3)/sqrt(sum(Y.y(:,3).^2));
YCE2 = Y.y(:,4)/sqrt(sum(Y.y(:,4).^2));
% 
% YS1 = Y.y(:,6)/sqrt(sum(Y.y(:,6).^2));
% YS2 = Y.y(:,7)/sqrt(sum(Y.y(:,7).^2));

CIn = CI/normCI;
CEn = CE/normCE;

figure
plot([YCI1 YCE1 CIn CEn])

YCI1'*CIn
YCE1'*CEn

YCI2'*CIn
YCE2'*CEn

figure
plot([YCI2 YCE2 CIn CEn])

%%

figure
plot(abs((fft(Y.y(:,1))./(fft(Y.y(:,7))+exp(-16)))))
hold on
plot(abs((fft(Y.y(:,3))./(fft(Y.y(:,7))+exp(-16)))),'r')

figure
plot(abs(fft(Ikernel)))



%%


M.opt.popFiringRates = 1;

[Ep,Cp,Eh,F,dFdp,dFdpp,f] = spm_nlsi_GN_LC_fr_spinit(M,U,Y)









