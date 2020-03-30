%% ML_FieldSimM
cd D:\mliete\MATLAB\FokkerPlanck_Special
% Simulation of cortical macroclumn with LIF neuron populations following
% Fokker-Planck dynamics.

%% initialization of pyramidal population

P.C  = .250; % nF => time seconds

P.Vl  = -57.5;  % mV
P.gl  = 250/20; % capacitance/time_constant [pF/ms] = [nS] => time seconds
P.Vr  = -70;    % mV 
P.Rt  = 3/1000; % s
P.d   = 0/1000; % Afferent axonal delay 
P.Vt  = -42;    % mV

P.G   = [1, 1];
P.Vg  = [-10,-70];
P.T   = [4,6]/1000; % REVIEW time constants
P.D   = 3000; % = 1/2 * sigma^2 ~ (mV^2)/ms
P.DL  = 0*[1, 1];
P.DGL = 0*[1, 3];

P.Vmin = -80;
P.Tr   = 40; % number of bins for Voltage space discretization
P.LQ   = 10; % number of bins for Queue space discretization
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

P.delay = min(P.Tr + floor(P.LQ*P.d/P.Rt),P.LVV-1);
P.FvarV   = FvarV;
P.Fsteady = Fsteady;
clear VV LVV T R FvarV Fsteady k m C

Ppyr = P;
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
P.T   = [4,6]/1000; % REVIEW time constants
P.D   = 3000; % = 1/2 * sigma^2 ~ (mV^2)/ms
P.DL  = 0*[1, 1];
P.DGL = 0*[1, 3];

P.Vmin = -80;
P.Tr   = 40; % number of bins for Voltage space discretization
P.LQ   = 10; % number of bins for Queue space discretization
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

P.delay = min(P.Tr + floor(P.LQ*P.d/P.Rt),P.LVV-1);
P.FvarV   = FvarV;
P.Fsteady = Fsteady;
clear VV LVV T R FvarV Fsteady k m C

Pii = P;

%% Initialization for Integration
np = 6;

StateSpace = pdf('norm',P.VV,-70,2);%pdf('norm',P.VV,-60,2);
%  StateSpace = ones(P.LVV,1);
StateSpace(P.Tr+1:end) = 0;
StateSpace = StateSpace/sum(StateSpace(:)); % probability density per milivolt
StateSpace = repmat(StateSpace,[1 1 np]);

N                = 10000;
SSoverT          = zeros(P.LVV,np,N);
SSoverT(:,:,1)   = StateSpace;
SSoverT(:,:,end) = StateSpace;
GV               = 0*ones(length(P.G),np,N);
Gout             = zeros(length(P.G),np);

% P.D   = 2800;
% P.DL  = 0*[1, 1]; %P.DL  = .5*[1 , 1];
% P.DGL = .00*[1/2, 1/2]; %P.DGL = .5*[1, 1];
% P.T   = 1*[4,16]/1000;

NodeP(1:np/2) = Ppyr;
NodeP(np/2+1:np) = Pii;
for l = 1:np
    SSoverT(:,l,[1 end]) = SSoverT(:,l,[1 end])/NodeP(l).Vres;
end

%% continue integration

dt = .0005;

% SSoverT(:,:,end) = SSoverT(:,:,k-1); GV(:,:,end) = GV(:,:,k-1);
SSoverT(:,:,1)   = SSoverT(:,:,end);
GV(:,:,1)        = GV(:,:,end);



U                = 2*ones(1,N) + .0*randn(1,N) + 0*.25*sin(4*2*pi*(dt:dt:N*dt)) + .0*(1/N:1/N:1);

np2=np/2;
GE = 1*[(1.0*(diag(ones(np2-1,1),1)+diag(ones(np2-1,1),-1))+.0*(diag(ones(np2-2,1),2)+diag(ones(np2-2,1),-2))) zeros(np2);
        2*eye(np2) zeros(np2)];
    
GI = 1*[zeros(np2) (diag(ones(np2,1)) + .0*(diag(ones(np2-1,1),1)+diag(ones(np2-1,1),-1))+.0*(diag(ones(np2-2,1),2)+diag(ones(np2-2,1),-2))) ;
        zeros(np2) zeros(np2)];
    
C = [1   zeros(1,np-1)]';

kini = 2;
%%
% figure(1)

% kini = k-1;
Drawplots = 0;
for k = kini:N       % integrate through time
   
    Gout(1,:) = GE*GV(1,:,k-1)' + C*U(k-1);
    Gout(2,:) = GI*GV(2,:,k-1)';
    
    for l = 1:np   % cycle through populations
        P = NodeP(l);
        SSS2  = [reshape(SSoverT(:,l,k-1),[],1,1); squeeze(GV(:,l,k-1))];
        [Qdpdt] = fx_LIFpopME(Gout(:,l),P);
        SSS            = expm(dt*Qdpdt)*SSS2;
        
        SSoverT(:,l,k) = SSS(1:P.LVV);
        GV(:,l,k)      = SSS(P.LVV+1:end);

%         subplot(np+1,1,l)
%         plot(P.VV,(SSoverT(:,l,k)));axis([P.Vmin P.Vmax -.001 .2]);
    end
    
    Gout(1,:) = .5*GE*(GV(1,:,k-1)'+ GV(1,:,k)') + .5*C*(U(k-1)+U(k));
    Gout(2,:) = .5*GI*(GV(2,:,k-1)'+ GV(2,:,k)');
    
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
% %     
    if Drawplots
        subplot(np+1,1,np+1)
        %     plot(dt:dt:N*dt,squeeze(GV(:,np,:))');
        subplot(np+1,1,np+1);plot(dt:dt:N*dt,squeeze(SSoverT(P.Tr+1,:,:)*P.Fsteady(P.Tr+1))')
        drawnow;
        disp(k*dt)
    else
        subplot(2,1,1)
        imagesc(SSoverT(:,:,k-1))
        subplot(2,1,2)
        plot(dt:dt:N*dt,squeeze(SSoverT(P.Tr+1,:,:)*P.Fsteady(P.Tr+1))')
        drawnow;
        disp(k*dt)
    end
%     input('')
    
end

% Plots
if Drawplots
    for l = 1:np
        P = NodeP(l);
        subplot(np+1,1,l);plot(P.VV,(SSoverT(:,l,k)));axis([P.VV(1) P.VV(end) -.001 .2]);drawnow;
    end
    subplot(np+1,1,np+1);plot(dt:dt:N*dt,squeeze(SSoverT(P.Tr+1,:,:)*P.Fsteady(P.Tr+1))')
else
    subplot(1,1,1)
    plot(dt:dt:N*dt,squeeze(SSoverT(P.Tr+1,:,:)*P.Fsteady(P.Tr+1))')
end
 %  subplot(np+1,1,np+1);plot(dt:dt:N*dt,squeeze(GV(:,np,:))')
    
      
%% LFP calculation: proportional to input currents for the pyramidal cell population

LFP = zeros(np,N);
GVI = zeros(size(GV)); % Input conductances from conductance outputs
GVI(1,:,:) = GE*squeeze(GV(1,:,:))+ C*U;
GVI(2,:,:)   = GI*squeeze(GV(2,:,:));

for l = 1:np
    P  =  NodeP(l);
    FV = P.FvarV;
    FV(P.Tr+1:end,:) = FV(repmat(P.R,P.LQ,1),:); % heuristic on refractory period
    LFP(l,:) = sum(squeeze(GVI(:,l,:)).*(-P.C*FV'*squeeze(SSoverT(:,l,:)*P.Vres))); % [nS*mV] = [pA] in pico Ampere per neuron (?)
end

figure
subplot(3,1,1);plot(dt:dt:N*dt,sum(LFP(1:np2,:))')
subplot(3,1,2);FFTplot(sum(LFP(1:np2,:))',1/dt)
subplot(3,1,3);TFanalysis(sum(LFP(1:np2,:)),7,100,1,250,1/dt,exp(1));

%%
% figure
pp=6;
subplot(3,1,1);plot((LFP(pp,:))')
subplot(3,1,2);FFTplot((LFP(pp,:))',1/dt)
subplot(3,1,3);TFanalysis((LFP(pp,:)),7,100,1,250,1/dt,exp(-4));

%%
SSoverTM = SSoverT/max(SSoverT(:));

for k = 2:N
    image(64*SSoverTM(:,:,k-1))
    drawnow
    disp(k)
end
%%
E = GE*squeeze(GV(1,:,:));
I = GI*squeeze(GV(2,:,:));

LFP2 = LFP(200:end);
%%

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
pp = 5;

OscillationTime = 557:2:577;% begining and ending of an oscillation cycle
Nisi = 300;
Nosc = length(OscillationTime);
resultsList = zeros(length(OscillationTime),Nisi);

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
            
            subplot(np+1,1,l)
            plot(P.VV,(SSoverT2(:,l,k)));axis([P.VV(1) P.VV(end) -.001 .2]);
        end
        
        subplot(np+1,1,np+1)
        subplot(np+1,1,np+1);plot(dt:dt:N*dt,squeeze(SSoverT2(P.Tr+1,:,:)*P.Fsteady(P.Tr+1))')
        drawnow;
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
ISIh(end)=0;
ISIh = ISIh/sum(ISIh)/dt;
bar(dt:dt:dt*300,ISIh')
hold on
plot(dt:dt:dt*Nisi,ISI,'r','LineWidth',2)
% hold off
xlim([dt,Nisi*dt])
ylim([0 .016/dt])
xlabel('\fontsize{12} Inter spike interval (s)')
ylabel('\fontsize{12} Probability density (s^-^1)')

