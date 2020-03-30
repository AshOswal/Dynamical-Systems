
cd D:\mliete\MATLAB\FokkerPlanck

P.C  = 10/1000; %42; % mF

P.Vl  = -70; % mV
P.gl  = 1; % mS
P.Vr  = -90; 
P.Rt  = 3/1000; %3; % s
P.d   = 2/1000;
P.Vt  = -40;

P.G   = [1, 1];
P.Vg  = [60,-90];
P.T   = [4,16]/1000;
P.D   = 3000; % = 1/2 * sigma^2 ~ (mV/ms)^2
P.DL  = 2*[1, 1];
P.DGL = 1*[1, 3];


P.Vmin = -100;
P.Vres = .1;
P.Vmax = P.Vt+10;

VV = P.Vmin:P.Vres:P.Vmax;
LVV= length(VV);

Fsteady = zeros(LVV,1);
FvarV   = zeros(LVV,length(P.G));
[C T]   = min(abs(VV-P.Vt)); 
[C R]   = min(abs(VV-P.Vr));

for k = 1:T
    Fsteady(k) = 1/P.C*(-P.gl*(VV(k)-P.Vl)) ;
    for m = 1:length(P.G)
        FvarV(k,m) = 1/P.C*(VV(k)-P.Vg(m));
    end
end

Fsteady(T+1:end) = (P.Vmax-P.Vt)/P.Rt;

P.VV    = VV;
P.LVV   = LVV;
P.Tr    = T;
P.delay = T + floor((LVV-T)/P.Rt*P.d);
P.R     = R;
P.FvarV   = FvarV;
P.Fsteady = Fsteady;
clear VV LVV T R FvarV Fsteady k m C

% Initialization for Integration
StateSpace = pdf('norm',P.VV,-50,3);%pdf('norm',P.VV,-54,8.3)
%  StateSpace = ones(P.LVV,1);
StateSpace(P.Tr+1:end) = 0;
StateSpace = (StateSpace/sum(StateSpace(:))/P.Vres); % probability density per milivolt
plot(P.VV,StateSpace)

G = zeros(size(P.G));

N=150;
SSoverT = zeros(P.LVV,N);
SSoverT(:,1) = StateSpace;
GV = 0*ones(length(P.G),N);%0*P.p/(1+P.p)*ones(length(P.G),N);
U  = randn(N,1); 

P.D   = 3000;
P.DL  = 10*[1, 1]; %P.DL  = .5*[1 , 1];
P.DGL = .01*[1/2, 1/2]; %P.DGL = .5*[1, 1];

dt = .001;

k=2;
StateSpace = [squeeze(SSoverT(:,k-1)); squeeze(GV(:,k-1))];
Go = [.4+.0*U(k) 0*GV(2,k-1)];

[Qdpdt] = fx_LIFpopME(Go,P);
%%


for k = 2:N
    
    StateSpace = [squeeze(SSoverT(:,k-1)); squeeze(GV(:,k-1))];
    Go = [.4+.0*U(k) 0*GV(2,k-1)];
    
    [Qdpdt] = fx_LIFpopME(Go,P);
    
    SSS          = spm_expm(dt*Qdpdt,StateSpace);
    SSoverT(:,k) = SSS(1:P.LVV);
    GV(:,k)      = SSS(P.LVV+1:end);
%     
    subplot(2,1,1)
    plot(P.VV,(SSoverT(:,k)),'k'); %alpha(.4);
    axis([P.Vmin P.Vmax -.0001 .2])
    xlabel('Voltage (mV)')
    ylabel('Probability Density (Neurons)')
    drawnow
    subplot(2,1,2)
    plot(GV')
    xlabel('Time (ms)')
    ylabel('Population Output (Spikes/Neuron/s)')
%     plot(P.VV,log(SSoverT(:,k)+exp(-16))); drawnow;
    disp([k sum(SSoverT(:,k))])
    
%     input('')
end

subplot(2,1,1)
plot(P.VV,(SSoverT(:,k))); %alpha(.4);
axis([P.Vmin P.Vmax -.0001 .2])
xlabel('Voltage (mV)')
ylabel('Probability Density (Neurons)')
drawnow
subplot(2,1,2)
plot(GV')
xlabel('Time (ms)')
ylabel('Population Output (Spikes/Neuron/s)')

%% Continue integration

SSoverT(:,1)=SSoverT(:,end);
SSoverT(:,2:end)=0;
GV(:,1) = GV(:,end);

%% Eigen values
% 
% P.D   = 1000;
% P.DL  = 0*[1, 2]; %P.DL  = .5*[1 , 1];
% P.DGL = 0*[1/2, 1/2]; %P.DGL = .5*[1, 1];
% Go = [.3 0];

[Qdpdt] = fx_LIFpopME(Go,P);

Qdpdt2 = Qdpdt(1:P.LVV,1:P.LVV);

[V S] = eig(full(Qdpdt2));
dS = diag(S);
[B,IX] = sort(real(diag(S)));
subplot(3,1,1)
plot(B)
subplot(3,1,2)
plot(real(dS),imag(dS),'.')
subplot(3,1,3)
plot(real([V(:,IX(end))]))
subplot(1,1,1)
plot([LIF_StSol(Go,P,0)' V(:,IX(end))/sum( V(:,IX(end)))*P.Vres])
%%
k=0
subplot(1,1,1)
plot(([real(V(:,IX(end:-1:end - k))) imag(V(:,IX(end:-1:end - k)))]))
% subplot(2,1,1)
% plot(([real(V(P.R:P.Tr,IX(end - k))) imag(V(P.R:P.Tr,IX(end - k)))]))
% subplot(2,1,2)
% A = fft(imag(V(P.R:P.Tr,IX(end - k))));
% plot(([abs(A(1:40)) real(A(1:40)) imag(A(1:40))]))
%%
StateSpace = pdf('norm',P.VV,-50,2);
StateSpace(P.Tr+1:end) = 0;
SS2 = StateSpace'/sum(StateSpace);
% %% 
% Saprox    = zeros(P.LVV,P.LVV);
% SaproxPhi = zeros(P.LVV,P.LVV);
% PhiMatrix = zeros(P.LVV,P.LVV);
% for k = 1:P.LVV/4
%    Saprox(:,k)    = V(:,IX(end:-1:end+1-k))*(V(:,IX(end:-1:end+1-k))\SS2);
%    aux            = cos(k*pi*(-P.LVV/2:P.LVV/2-1)/P.LVV)'+1i*sin(k*pi*(-P.LVV/2:P.LVV/2-1)/P.LVV)';
%    PhiMatrix(:,k) = V(:,IX(end)).*aux/norm(aux);
%    SaproxPhi(:,k) = PhiMatrix(:,end:-1:end+1-k)*(PhiMatrix(:,end:-1:end+1-k)\SS2);
%    disp([k P.LVV])
% end
% %%
% figure
% subplot(2,1,1)
% plot(real(Saprox(:,1:end)))
% subplot(2,1,2)
% plot(real(SaproxPhi(:,1:end)))



%%
Nb = 20;

Nb = 10;
x = pi*((1:P.LVV)'-P.R)/(P.LVV-P.R);
Bsin = zeros(P.LVV,Nb);
Bsin2 = zeros(P.LVV,Nb);
for k = 1:Nb
    Bsin(:,k) = (cos(k*x) + 1i*sin(k*x));
    Bsin(:,k)  = Bsin(:,k)/norm(Bsin(:,k));
    Bsin2(:,k) = (x.^(k));
    Bsin2(:,k) = Bsin2(:,k)/norm(Bsin2(:,k)); 
end
Bsin = [ones(P.LVV,1) Bsin];
Bsin2 = [ones(P.LVV,1) Bsin2];
Nb = Nb+1;

% plot(real(Bsin))
%%

QB = Bsin\Qdpdt2*Bsin;
[Vb Sb] = eig(full(QB));
dSb = diag(Sb);
[Bb,IXb] = sort(real(diag(Sb)));
subplot(3,1,1)
plot(Bb)
subplot(3,1,2)
plot(real(dSb),imag(dSb),'.')
subplot(3,1,3)
plot(real([Vb(:,IXb(end))]))

%%
N = 1000;
RandBF = zeros(P.LVV,N);
flag = .73 ;

for k = 1:N
    Go = [2*rand; 2*rand];
    RandBF(:,k) = LIF_StSol(Go,P,flag);
    disp(k)
end

plot(RandBF)
RandBF(isnan(RandBF)) = 0;
RandBF(isinf(RandBF)) = 0;
[U,S,~] = svd(RandBF,0);
plot(U(:,1:Nb)*(U(:,1:Nb)\RandBF))
BFU = U(:,1:Nb);
plot(BFU)
%% Integration using the different basis functions

Bsin2 = BFU;

N = 500;
QB = Bsin\(Qdpdt2*Bsin);
[Vb Sb] = eig(full(QB));

QB2 = Bsin2\(Qdpdt2*Bsin2);
[Vb2 Sb2] = eig(full(QB2));

SSoverT  =zeros(701,N);
SSoverT(:,1) = StateSpace;

piB  = zeros(Nb,N);
piB2 = zeros(Nb,N);
piE  = zeros(701,N);

piB(:,1)  = Vb\(Bsin\SSoverT(:,1));
piB2(:,1) = Vb2\(Bsin2\SSoverT(:,1));
piE(:,1)  = V\SSoverT(:,1);

dt = .001;

for k = 2:N;
    piB(:,k)  = exp(k*dt*diag(Sb)).*piB(:,1);
    piB2(:,k) = exp(k*dt*diag(Sb2)).*piB2(:,1);
    piE(:,k)  = exp(k*dt*diag(S)).*piE(:,1);
end


SSoverTB  = Bsin*Vb*piB;
SSoverTB2 = Bsin2*Vb2*piB2;
SSoverT   = V*piE;

%%
subplot(1,1,1)
k = 1
plot([real(SSoverT(:,k)) real(SSoverTB(:,k)) real(SSoverTB2(:,k)) ])
%  ylim([-0.05 .25])

%%
 plot([real(SSoverT(end,:)') real(SSoverTB(end,:))' real(SSoverTB2(end,:))' ])
%%
subplot(1,1,1)
k = 73
fSS=real(fft(real(SSoverT(:,k))));
plot(fSS(1:50))
 ylim([-4 10])

%%
k = -8;
TestA=.001*(Qdpdt2*Bsin2);
plot([abs(TestA(:,k)) abs(Bsin(:,k))])




