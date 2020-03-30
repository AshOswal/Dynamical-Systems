%% ML_NodeSimM
clear all
cd /Users/mfpleite/Documents/MATLAB/FokkerPlanck_Special
% Simulation of cortical macroclumn with LIF neuron populations following
% Fokker-Planck dynamics.

%% initialization of pyramidal population

P.C  = .250; % [nF] => time seconds

P.Vl  = -57.5;  % [mV]
P.gl  = 250/20; % capacitance/time_constant [pF/ms] = [nS] => time seconds
P.Vr  = -70;    % [mV] reset voltage
P.Rt  = 3/1000; % [s]  refractory time
P.d   = 0/1000; % [s]  afferent axonal delay 
P.Vt  = -42;    % [mV] spike threshold voltage
P.G   = [1, 1]; % [1]  synaptic eficacy scaling
P.Vg  = [-10,-70]; % [mV] reversal portential of the families of synaptic channels 
P.T   = [4,16]/1000; % [s] decay time constants for the synaptic conductances REVIEW
P.D   = 6000; % = 1/2 * sigma^2 ~ [(mV^2)/ms] membrane noise term
P.DL  = 0*[1, 1]; % Unused parameters to model multiplicative membrane noise
P.DGL = 0*[1, 3]; % Unused parameters to model multiplicative membrane noise

P.Vmin = -80; % [mV] minimum of the voltage space
P.Tr   = 80; % number of bins for Voltage space discretization
P.LQ   = 20; % number of bins for Queue space discretization


P.LVV  = P.Tr+P.LQ;
P.VV   = linspace(P.Vmin,P.Vt,P.Tr); %P.Vt*ones(1,P.LQ)];
P.Vres = P.VV(2)-P.VV(1);
P.VV   = [P.VV (P.VV(end)+P.Vres):P.Vres:(P.VV(end)+P.LQ*P.Vres)];

Fsteady = zeros(P.LVV,1);
FvarV   = zeros(P.LVV,length(P.G));
[C R]   = min(abs(P.VV-P.Vr));
P.R     = R; % reset bin - a small discretization error is made here

for k = 1:P.Tr
    Fsteady(k) = 1/P.C*(-P.gl*(P.VV(k)-P.Vl)) ;
    for m = 1:length(P.G)
        FvarV(k,m) = 1/P.C*(P.VV(k)-P.Vg(m));
    end
end

Fsteady(P.Tr+1:end) = P.LQ/P.Rt;

P.delay = min(P.Tr + 1 + floor(P.LQ*P.d/P.Rt),P.LVV-1);
P.FvarV   = FvarV;
P.Fsteady = Fsteady;
clear VV LVV T R FvarV Fsteady k m C

Ppyr = P;

% [Qdpdt] = fx_LIFpopME([20;10],P);
% 
% [V S] = eig(full(Qdpdt(1:end-2,1:end-2)));
% [B,IX] = sort(real(diag(S)));
% plot(real(diag(S)),imag(diag(S)),'.')

%% initialization of inhibitory interneurons

P.C  = .095; % nF => time seconds

P.Vl  = -65;    % mV
P.gl  = 95/10;  % capacitance/time_constant [pF/ms] = [nS] => time seconds
P.Vr  = -70;    % mV 
P.Rt  = 3/1000; % s
P.d   = 0/1000; % Afferent axonal delay 
P.Vt  = -45;    % mV

P.G   = [1, 1];
P.Vg  = [-10,-70];
P.T   = [4,16]/1000; % REVIEW time constants
P.D   = 6000; % = 1/2 * sigma^2 ~ (mV^2)/ms
P.DL  = 0*[1, 1];
P.DGL = 0*[1, 3];

P.Vmin = -80;
P.Tr   = 80; % number of bins for Voltage space discretization
P.LQ   = 20; % number of bins for Queue space discretization
P.LVV  = P.Tr+P.LQ;
P.VV   = linspace(P.Vmin,P.Vt,P.Tr); %P.Vt*ones(1,P.LQ)];
P.Vres = P.VV(2)-P.VV(1);
P.VV   = [P.VV (P.VV(end)+P.Vres):P.Vres:(P.VV(end)+P.LQ*P.Vres)];

Fsteady = zeros(P.LVV,1);
FvarV   = zeros(P.LVV,length(P.G));
[C R]   = min(abs(P.VV-P.Vr));
P.R     = R; % reset bin - a small discretization error is made here

for k = 1:P.Tr
    Fsteady(k) = 1/P.C*(-P.gl*(P.VV(k)-P.Vl)) ;
    for m = 1:length(P.G)
        FvarV(k,m) = 1/P.C*(P.VV(k)-P.Vg(m));
    end
end

Fsteady(P.Tr+1:end) = P.LQ/P.Rt;

P.delay = min(P.Tr + 1 + floor(P.LQ*P.d/P.Rt),P.LVV-1);
P.FvarV   = FvarV;
P.Fsteady = Fsteady;
clear VV LVV T R FvarV Fsteady k m C

Pii = P;

%% Initialization for Integration
np = 2;

StateSpace = pdf('norm',P.VV,-70,2);%pdf('norm',P.VV,-60,2);
%  StateSpace = ones(P.LVV,1);
StateSpace(P.Tr:end) = 0;
StateSpace = StateSpace/sum(StateSpace(:)); % probability density per milivolt
StateSpace = repmat(StateSpace,[1 1 np]);

N                = 5000;
SSoverT          = zeros(P.LVV,np,N);
SSoverT(:,:,1)   = StateSpace;
SSoverT(:,:,end) = StateSpace;
GV               = 0*ones(length(P.G),np,N);
Gout             = zeros(length(P.G),np);

if np == 2;
    NodeP(1) = Pii;
    NodeP(2) = Ppyr;
elseif np ==3
    NodeP(1) = Pii;
    NodeP(2) = Pii;
    NodeP(3) = Ppyr;
elseif np==4
    NodeP(1) = Pii;
    NodeP(2) = Ppyr;
    NodeP(3) = Pii;
    NodeP(4) = Ppyr;
end
for l = 1:np
    SSoverT(l,:,[1 end]) = SSoverT(l,:,[1 end])/NodeP(l).Vres;
end

LFP = zeros(np,N);
GVI = zeros(size(GV)); % Input conductances from conductance outputs


kini = 2;
%% continue integration

dt = .0005;

SSoverT(:,:,1)   = SSoverT(:,:,end);
GV(:,:,1)        = GV(:,:,end);
GVI(:,:,1)        = GVI(:,:,end);
U                = [2.5*ones(1,N); 0*ones(1,N);];

if np ==3
    GE = 1*[0    0   .25;
            0    0   1;
            1    0   0];
    
    GI = 1*[0   .5  0;
            0   0   0;
            0   1   0];
    
    C = [ 1 0 0]';
elseif np==2
    GE = 3*[ 0    4;
             0    .3];
    
    GI = .5*[.125   0;
            .75   0];
    
    C = [.5 1]';
elseif np==4
    GE = .5*[ 0    1    0   0;
             0    .5    0   0;
             0     0    0   1;
             0    .0    0   0,];
    
    GI = 1*[.2   0   0   0;
             1   0   0   0;
             0   0   .2   0;
             0   0   1   0];
    
    C = [0 0 0 1]';

end
% 
% GE = .030*[0    0   1;
%           0    0   .25;
%          .8   0    0  ];
% 
% GI = .0475*[0  0*1/4  0;
%           0   0   0;
%           0   1   0];
% C = [ 0 0 1]';

% GE   =   1*[ 0     0     0     0
%              1     0     0     0
%              1     0     0     1
%              0     1     0     0];
%  
% GI   =   1*[ 0     0     1     0
%              0     0     0     0
%              0     0     0     0
%              0     1     1     0];
%  C = [ 1 0 0 0]';
if kini==N kini=2; end

%%
% figure(1)
Drawplots = 1;
RTcontrol = zeros(size(U,1),1);
for k = kini:N       % integrate through time    
    kini = k;
%     input('press enter')
    % real time control
    RTtc = 6/1000;
    RTgain = 0*10^-4;
    RTcontrol = RTcontrol + dt/RTtc*(RTgain*[max((LFP(2,k-1)+0*LFP(2,max(k-2,1)))/dt,0);0] - RTcontrol);
    U(:,k-1)  = U(:,k-1) + RTcontrol;
    
    % --- Prediction step ---
    Gout(1,:) = GE*GV(1,:,k-1)';
    Gout(2,:) = GI*GV(2,:,k-1)';
    Gout      = Gout + U(:,k-1)*C';
    
    for l = 1:np   % cycle through populations
        P = NodeP(l);
        SSS2  = [reshape(SSoverT(:,l,k-1),[],1,1); squeeze(GV(:,l,k-1))];
        [Qdpdt] = fx_LIFpopME(Gout(:,l),P);
        SSS            = expm(dt*Qdpdt)*SSS2;
        
        SSoverT(:,l,k) = SSS(1:P.LVV);
        GV(:,l,k)      = SSS(P.LVV+1:end);
    end
    
    % --- Correction step ---
    Gout(1,:) = .5*GE*(GV(1,:,k-1)'+ GV(1,:,k)');
    Gout(2,:) = .5*GI*(GV(2,:,k-1)'+ GV(2,:,k)');
    Gout      = Gout + .5*(U(:,k-1)+U(:,k))*C';
    
    for l = 1:np   % cycle through populations
        P = NodeP(l);
        SSS2    = [reshape(SSoverT(:,l,k-1),[],1,1); squeeze(GV(:,l,k-1))];
        [Qdpdt] = fx_LIFpopME(Gout(:,l),P);
        SSS     = expm(dt*Qdpdt)*SSS2;
        
        SSoverT(:,l,k) = SSS(1:P.LVV);
        GV(:,l,k)      = SSS(P.LVV+1:end);
        if Drawplots
            subplot(np+1,1,l)
            plot(P.VV,(SSoverT(:,l,k)));axis([P.VV(1) P.VV(end) -.001 .2]);
            drawnow;
        end
    end
    
    % --- LFP computation ---   
    GVI(:,:,k)= Gout; %- .5*(U(:,k-1)+U(:,k))*C';
    for l = 1:np
        P  =  NodeP(l);
        FV = P.FvarV;
        FV(P.Tr+1:end,:) = FV(repmat(P.R,P.LQ,1),:); % heuristic on refractory period ignoring repolarization currents
        LFP(l,k) = sum(GVI(:,l,k).*(P.C*FV'*squeeze(SSoverT(:,l,k)*P.Vres))); % [nS*mV] = [pA] in pico Ampere per neuron (?)
    end

    
    
    if Drawplots
        subplot(np+1,1,np+1)
        plot(dt:dt:N*dt,LFP(2,:));
%         plot(dt:dt:N*dt,70*squeeze(GVI(:,2,:))');
%         G=zeros([size(GE) 2]);G(:,:,1) = GE; G(:,:,2) = GI; U2 =  U(:,1)*C'; S2 = [SSoverT(:,:,kini-1); squeeze(GV(:,:,kini-1))];
%         [~, J] = fx_LIFpopMEJ(NodeP,G,U2,S2);
%         [~, S1] = eig(full(J));
%         plot(real(diag(S1)),imag(diag(S1)),'.')
%         xlim([-5 1]*10^4); ylim([-1 1]*10^4)
%         subplot(np+1,1,np+1);plot(dt:dt:N*dt,squeeze(SSoverT(P.Tr+1,:,:)*P.Fsteady(P.Tr+1))')
        drawnow;
        disp(k*dt)
    else
        disp(['Timestep: ' num2str(k*dt*1000) ' ms'])
    end   
end

% Plots
 for l = 1:np
     P = NodeP(l);
     subplot(np+1,1,l);plot(P.VV,(SSoverT(:,l,k)));axis([P.VV(1) P.VV(end) -.001 .2]);drawnow;
 end
 subplot(np+1,1,np+1);plot(dt:dt:N*dt,squeeze(SSoverT(P.Tr+1,:,:)*P.Fsteady(P.Tr+1))')
%  subplot(np+1,1,np+1);plot(dt:dt:N*dt,squeeze(GV(:,np,:))')
figure;plot(squeeze(GVI(1,2,:))*70,squeeze(GVI(2,2,:))*70)    
%%

G(:,:,1) = GE; G(:,:,2) = GI; U2 =  U(:,1)*C'; S2 = [SSoverT(:,:,kini-1); squeeze(GV(:,:,kini-1))];
S0=S2;
iterN = 0;
f = Inf;
while norm(f)>.001
    iterN = iterN +1; disp(iterN)
    [f J] = fx_LIFpopMEJ(NodeP,G,U2,S0);
    % constrain x to sum 1! otherwise J is singular
    f = [f;zeros(np,1)]; 
    J = [J;ones(1,100) zeros(1,104);zeros(1,102) ones(1,100) 0 0];
%     J = [J;ones(1,100) zeros(1,206);zeros(1,102) ones(1,100) 0 0 zeros(1,102);
%         zeros(1,204) ones(1,100) 0 0];
    S0 = S0(:) - J\f;
    S0 = reshape(S0,size(S2));
    plot(S0)
end
SSoverT(:,:,kini-1) = S0(1:100,:); GV(:,:,kini-1)      = S0(101:end,:);


%% LFP calculation: proportional to input currents for the pyramidal cell population

LFP = zeros(np,N);
GVI = zeros(size(GV)); % Input conductances from conductance outputs
GVI(1,:,:) = GE*squeeze(GV(1,:,:))+ C*U(1,:);
GVI(2,:,:)   = GI*squeeze(GV(2,:,:))+ C*U(2,:);

for l = 1:np
    P  =  NodeP(l);
    FV = P.FvarV;
    FV(P.Tr+1:end,:) = FV(repmat(P.R,P.LQ,1),:); % heuristic on refractory period
    LFP(l,:) = sum(squeeze(GVI(:,l,:)).*(P.C*FV'*squeeze(SSoverT(:,l,:)*P.Vres))); % [nS*mV] = [pA] in pico Ampere per neuron (?)
end

figure
subplot(3,1,1);plot(LFP')
subplot(3,1,2);FFTplot(LFP(l,:)',1/dt)
subplot(3,1,3);TFanalysis(LFP(l,:),7,100,1,450,1/dt,exp(0));



%%
P = NodeP(2);
[Qdpdt] = fx_LIFpopME([7;10],P);

[V S] = eig(full(Qdpdt(1:end-2,1:end-2)));
[VT ST] = eig(full(Qdpdt(1:end-2,1:end-2)'));
[B,IX] = sort(real(diag(S)));
[BT,IXT] = sort(real(diag(ST)));
B = diag(S(IX,IX));
BT = diag(ST(IXT,IXT));
plot(real(V(:,IX(end))))

VS = V(:,IX);
VST = VT(:,IXT);

%% ISI

pp = 2;

OscillationTime = 157:2:177;% begining and ending of an oscillation cycle
Nisi = 300;
Nosc = length(OscillationTime);
resultsList = zeros(length(OscillationTime),Nisi);

Drawplots = 1;
for m = 1:Nosc
    initialT = OscillationTime(m);
    SSoverT2 = zeros(size(SSoverT));
    SSoverT2(P.R,:,initialT-1) = 1/P.Vres;
    for k = initialT:initialT+Nisi       % integrate through time
        
        Gout(1,:) = GE*GV(1,:,k-1)' + C*U(k-1);
        Gout(2,:) = GI*GV(2,:,k-1)';
        
        for l = pp   % cycle through populations
            
            SSS2  = [reshape(SSoverT2(:,l,k-1),[],1,1); squeeze(GV(:,l,k-1))];
            [Qdpdt] = fx_LIFpopME_ISI(Gout(:,l),P);
            SSS            = expm(dt*Qdpdt)*SSS2;            
            SSoverT2(:,l,k) = SSS(1:P.LVV);
            
            if Drawplots
                subplot(np+1,1,l)
                plot(P.VV,(SSoverT2(:,l,k)));axis([P.VV(1) P.VV(end) -.001 .2]);
            end
        end
        
        if Drawplots
            subplot(np+1,1,np+1)
            subplot(np+1,1,np+1);plot(dt:dt:N*dt,squeeze(SSoverT2(P.Tr+1,:,:)*P.Fsteady(P.Tr+1))')
            drawnow;
        end
    end
    
    resultsList(m,:) = SSoverT2(P.Tr+1,pp,OscillationTime(m):(OscillationTime(m)+Nisi-1));
    disp(['Time: ' num2str(m) ' out of ' num2str(Nosc)])
end

subplot(np+1,1,np+1)
hold on
plot(dt:dt:N*dt,squeeze(SSoverT(P.Tr+1,:,:)*P.Fsteady(P.Tr+1))','k')
hold off

%%
figure
ISI = resultsList'*squeeze(SSoverT(P.Tr+1,pp,(OscillationTime)-P.Rt/dt-1));
ISI = ISI/sum(ISI)/dt;
try
ISIh(end)=0;
ISIh = ISIh/sum(ISIh)/dt;
bar(dt:dt:dt*300,ISIh')
hold on
end
plot(dt:dt:dt*Nisi,ISI,'r','LineWidth',2)
% hold off
xlim([dt,Nisi*dt])
ylim([0 .016/dt])
xlabel('\fontsize{12} Inter spike interval (s)')
ylabel('\fontsize{12} Probability density (s^-^1)')


