%% ML_NodeSimM
clear all
cd '/home/mleite/Documents/MATLAB/FokkerPlanck_Special'
% Simulation of cortical macroclumn with LIF neuron populations following
% Fokker-Planck dynamics.
%%
% P = NodeP(2);
[Qdpdt] = fx_LIFpopME([0;0],Ppyr);

[V S] = eig(full(Qdpdt(1:end,1:end)));
[VT ST] = eig(full(Qdpdt(1:end,1:end)'));
[B,IX] = sort(real(diag(S)));
[BT,IXT] = sort(real(diag(ST)));
B = diag(S(IX,IX));
BT = diag(ST(IXT,IXT));
pdfs = real(V(:,IX(end)));
pdfs = pdfs/sum(pdfs(1:end-2))/P.Vres;
figure;plot(pdfs(1:end-2))
disp(pdfs(end))
pdflist = {pdflist pdfs};
VS = V(:,IX);
VST = VT(:,IXT);

figure;
plot(squeeze(GVI(1,2,:))*70,squeeze(GVI(2,2,:))*70);
xlabel('GE');ylabel('GI')
hold on
plot(squeeze(GVI2(1,2,:))*70,squeeze(GVI2(2,2,:))*70,'k');
xlabel('GE');ylabel('GI')
plot(squeeze(GVI1(1,2,:))*70,squeeze(GVI1(2,2,:))*70,'r');
xlabel('GE');ylabel('GI')




%% initialization of pyramidal population ---------------------------------
% parameters for pyramidal neurons from:
% Badel,L. et al.(2008) BiolCybern 99,361-370
% Badel,L. et al.(2008) J Neurophysiol 99, 656-666
% -------------------------------------------------------------------------
r = 1;
P.C   = .250;           % [nF] Cell capacitance
P.Vl  = -50.5;          % [mV] Leak reversal portential
P.gl  = P.C/0.02;       % [pF/ms] = [nS] Leak conductance�
P.Rt  = 3/1000;         % [s]  Refractory time
P.Vg  = [0,-70];        % [mV] Reversal portential of the families of synaptic channels 
P.T   = [10,16]/1000;   % [s] Decay time constants for the synaptic conductances
P.D   = 6000;           % (1/2*sigma^2) [(mV^2)/s] Membrane noise term

P.Vt  = -42;            % [mV] Spike threshold voltage
P.Vr  = -70;            % [mV] Reset voltage
P.d   = 0/1000;         % [s]  Afferent axonal delay 
P.G   = [1, 1];         % [1]  Synaptic eficacy scaling (not in use)
P.DL  = 0*[1, 1];       % Unused parameters to model multiplicative membrane noise
P.DGL = 0*[1, 3];       % Unused parameters to model multiplicative membrane noise

P.Vmin = -80;           % [mV] minimum of the voltage space
P.Tr   = r*80;          % number of bins for Voltage space discretization
P.LQ   = r*20;          % number of bins for Queue space discretization

% -------- Do not Change --------
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
% -------- \Do not Change --------

Ppyr = P; % save final structure

%% initialization of inhibitory interneurons ------------------------------
% parameters for inhibitory interneuons from:
% Badel,L. et al.(2008) BiolCybern 99,361-370
% Badel,L. et al.(2008) J Neurophysiol 99, 656-666
% -------------------------------------------------------------------------

P.C   = .095;           % [nF] Cell capacitance
P.Vl  = -65;            % [mV] Leak reversal portential
P.gl  = 95/10;          % [pF/ms] = [nS] Leak conductance
P.Rt  = 3/1000;         % [s]  Refractory time
P.Vg  = [0,-70];        % [mV] Reversal portential of the families of synaptic channels
P.T   = [10,16]/1000;   % [s] Decay time constants for the synaptic conductances
P.D   = 6000;           % (1/2*sigma^2) [(mV^2)/ms] Membrane noise term

P.Vt  = -45;            % [mV] Spike threshold voltageCtrl+X, Ctrl+S
P.Vr  = -70;            % [mV] Reset voltage
P.d   = 0/1000;         % [s]  Afferent axonal delay 
P.G   = [1, 1];         % [1]  Synaptic eficacy scaling (not in use)
P.DL  = 0*[1, 1];       % Unused parameters to model multiplicative membrane noise
P.DGL = 0*[1, 3];       % Unused parameters to model multiplicative membrane noise

P.Vmin = -80;           % [mV] minimum of the voltage space
P.Tr   = r*80;          % number of bins for Voltage space discretization
P.LQ   = r*20;          % number of bins for Queue space discretization

% -------- Do not Change --------
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
% -------- \Do not Change --------

Pii = P; % save final structure

%% Initialization for Integration -----------------------------------------
% run to reset all simultation related variables
% -------------------------------------------------------------------------

np = 2; % number of populations

StateSpace = pdf('norm',P.VV,-70,2);        % initializes populations with gaussian probability density functions (pdfs)
StateSpace(P.Tr:end) = 0;                   % clears queue
StateSpace = StateSpace/sum(StateSpace(:)); % normalization to probability density per milivolt
StateSpace = StateSpace/P.Vres;
StateSpace = repmat(StateSpace,[1 1 np]);   % All populations start from the same distributions

N                = 1000;                    % Number of time bins to integrate
SSoverT          = zeros(P.LVV,np,N);       % Population pdfs over time
SSoverT(:,:,1)   = StateSpace;              % Initialization
SSoverT(:,:,end) = StateSpace;

GV               = 0*ones(length(P.G),np,N); % Afferent conductances over time
Gout             = zeros(length(P.G),np);    

% ----- Definition of the structures for each population -----
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
    SSoverT(l,:,[1 end]) = SSoverT(l,:,[1 end])/NodeP(l).Vres; % pdf normalization for Vres
end
% -------------------------------------------------------------

LFP = zeros(np,N);
GVI = zeros(size(GV));              % Input conductances from conductance outputs
U = [1*ones(1,N); 0*ones(1,N);];    % external input conductance to each channel family
RTcontrol = ones(1,size(U,2));      % Real time control from LFP values (under development)

kini = 2;                           % reset integration counter

%% continue integration from the previous position ------------------------
% -------------------------------------------------------------------------

dt = .0001;         % time step (for most purposes .0005 seconds is more than enough)

SSoverT(:,:,1)   = SSoverT(:,:,end); % continue from last iteration
GV(:,:,1)        = GV(:,:,end);
GVI(:,:,1)       = GVI(:,:,end);


U = [1*ones(1,N); 0*ones(1,N);]; % external input conductance to each channel family

% ---- Example connectivity matrices between populations ------------------
if np ==3
    GE = 1*[0    0   .25;
            0    0   1;
            1    0   0];
    
    GI = 1*[0   .5  0;
            0   0   0;
            0   1   0];
    
    C = [ 1 0 0]';
elseif np==2 % CA3 example
    GE = 1*[ 0    1;
             0    .3];
    
    GI = 1*[.125   0;
            .75   0];
    
    C = [.0 3]';
elseif np==4 % CA3 CA1 example
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
% -------------------------------------------------------------------------

if kini==N kini=2; end
RTcontrol(1:2) = RTcontrol(end-1:end); % Real time control from LFP values (under development)

%% Integrate over time ----------------------------------------------------
% The integration can be interrupted using ctrl+c and then continue from
% the same place, with or without changes to the parameters.
% -------------------------------------------------------------------------

Realtimecontrol = 0;
Drawplots = 1;                  % set to 1 or 0 to dysplay graphics or not
LocalElectrode = 1;

for k = kini:N % integrate through time    
    kini = k;
    
    %--- Real time control - under development
    if Realtimecontrol
        RTtc = 5/1000;
        RTgain =  00*10^-4;
        RTgaindt = 0;
        RTcontrol(1,k) = RTcontrol(1,k-1) + dt/RTtc*(...
            max(1 + RTgain*(LFP(2,k-1) - LFP(2,max(k-2,1))),0)...
            - RTcontrol(1,k-1));
        U(:,k-1)  =  RTcontrol(1,k);
        U(:,k)    =  RTcontrol(1,k);
    end
    % ---------------------------------------------------
    
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
            subplot(np+2,1,l)
            plot(P.VV,(SSoverT(:,l,k)));axis([P.VV(1) P.VV(end) -.001 .2]);
            drawnow;
        end
    end
    
    % --- LFP computation ---   
    GVI(:,:,k)= Gout; % Synaptic input conductances received;
    for l = 1:np
        P  =  NodeP(l);
        FV = P.FvarV;
        FV(P.Tr+1:end,:) = FV(repmat(P.R,P.LQ,1),:); % heuristic on refractory period
        LFP(l,k) = [1 1]*(GVI(:,l,k).*(FV'*squeeze(SSoverT(:,l,k)*P.Vres))); % [nS*mV] = [pA] in pico Ampere (per neuron)
        if LocalElectrode
            LFP(l,k) = LFP(l,k) + squeeze(SSoverT(P.Tr+1,l,k)*P.Vres)*P.C*(P.Vt-P.Vr)/.001; % Current for repolarization in 1 ms (per neuron)
        end
    end
    
    %     --- Plots ----
    
    pop = 2; % Pyramidal cells
    P   =  NodeP(pop);
    EIcurrents = [squeeze(GVI(1,pop,:))*(P.Vg(1)-P.Vg(2)) + P.gl*(P.Vl-P.Vg(2)),...
        squeeze(GVI(2,pop,:))*(P.Vg(1)-P.Vg(2)) + P.gl*(P.Vg(1)-P.Vl)];
    
    if Drawplots
        subplot(np+2,1,np+1)
        plot(dt:dt:N*dt,LFP(pop,:),'k');
        hold on
        plot(dt:dt:N*dt,EIcurrents)
        hold off
        
        subplot(np+2,1,np+2)
        plot(EIcurrents(:,1),EIcurrents(:,2));
        xlabel('E (pA)'); %xlim([-100 100]+EIcurrents(k,1)); 
        ylabel('I (pA)'); %ylim([-100 100]+EIcurrents(k,2));
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
figure;plot(squeeze(GVI(1,2,:))*70,squeeze(GVI(2,2,:))*70);xlabel('GE');ylabel('GI')

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

figure
subplot(3,1,1);plot(LFP')
subplot(3,1,2);FFTplot(LFP(l,:)',1/dt)
subplot(3,1,3);TFanalysis(LFP(l,:),7,100,1,450,1/dt,exp(0));



%%
P = NodeP(2);
[Qdpdt] = fx_LIFpopME([0;0],Ppyr);

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

OscillationTime = 217:2:277;% begining and ending of an oscillation cycle
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
xlabel('\fontsize{24} Inter spike interval (s)')
ylabel('\fontsize{24} Probability density (s^-^1)')


