function [ dxdp, x ] = spm_lifpopsys_LC_adapstep(P,M,Vp)
%   spm_lifpopsys_LC integrates a system of LIF populations
% 
%   P - structure of system parameters with the following fields:
% 
%       P.A - Connectivity matrices for each synaptic channel family;
%       P.C - External input matrices for each synaptic channel family;
%       P.P - Variable population parameters;
%
%   M - structure defining numerical parameters and integration options:
%
%       M.P   - Parameter structures for each population extended with fixed
%               parameters and auxiliary fields for numerical integration;
%       M.opt - Options for numerical integration;
%
%   Vp - Matrix of directions of P w.r.t. wich differentiate x.
%
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging



%% Test Parameters
%
clear all

P.ang    = -pi/2;

P.A(1).M = [ -32  log(1);
             -32  log(.5)];
    
P.A(2).M = [log(.5) -32;
            log(1) -32];
    
P.C(1).M = [-32 log(3)]';

P.C(2).M = [-32 -32]';

% --- inhibitory interneurons ---
P.P(1).C   = log(.095);           % [nF] Cell capacitance
P.P(1).Vl  =    -65;              % [mV] Leak reversal portential
P.P(1).gl  = log(95/10);          % [pF/ms] = [nS] Leak conductance
P.P(1).Rt  = log(3/1000);         % [s]  Refractory time
P.P(1).Vg  = [0,-70];             % [mV] Reversal portential of the families of synaptic channels
P.P(1).T   = log([10,16]/1000);   % [s] Decay time constants for the synaptic conductances
P.P(1).D   = log(6000);           % (1/2*sigma^2) [(mV^2)/ms] Membrane noise term

% --- excitatory pyramidal cells ---
P.P(2).C   = log(.250);            % [nF] Cell capacitance
P.P(2).Vl  =    -50.5;             % [mV] Leak reversal portential
P.P(2).gl  = P.P(2).C - log(0.02); % [pF/ms] = [nS] Leak conductance�
P.P(2).Rt  = log(3/1000);          % [s]  Refractory time
P.P(2).Vg  = [0,-70];              % [mV] Reversal portential of the families of synaptic channels 
P.P(2).T   = log([10,16]/1000);    % [s] Decay time constants for the synaptic conductances
P.P(2).D   = log(6000);            % (1/2*sigma^2) [(mV^2)/s] Membrane noise term

M.opt.N = 1000;
M.r     = 1;

Vp = 0*spm_vec(P); Vp = [Vp Vp]; Vp(4,1) = 1; Vp(4,2) = -1;

% Initialization for Integration -----------------------------------------
% -------------------------------------------------------------------------
 
if ~spm_lifpopsys_LC_check(P,M,Vp); return; end; % che+ck inputs

P0  = P;                % Save original parameters
P   = spm_ExpP(P0);
M   = spm_lifpopsys_LC_prepare(P,M);
M.opt.Drawplots = 1;



%% Integrate over time ----------------------------------------------------
% The integration can be interrupted using ctrl+c and then continue from
% the same place, with or without changes to the parameters.
% -------------------------------------------------------------------------

Gout       = zeros(M.nc,M.np);
M.LFP      = zeros(M.np,M.opt.N);
M.Currents = zeros(M.np,M.nc,M.opt.N);


if ~isfield(M.opt,'Drawplots')
    M.opt.Drawplots = 0; 
elseif M.opt.Drawplots
    figure(M.opt.Drawplots)
end


if ~isfield(M.opt,'dttol'); M.opt.dttol = 10^-3; end;
if ~isfield(M.opt,'dt'); M.opt.dt = .001; end;
if ~isfield(M.opt,'T0'); M.opt.T0 = .12; end; 
if ~isfield(M.opt,'svdtol'); M.opt.svdtol = 10^-6; end; 
if ~isfield(M.opt,'dpsize'); M.opt.dpsize = exp(-3); end;
if ~isfield(M.opt,'Ncoefs'); M.opt.Ncoefs = 1000; end;

M.Gv0 = zeros(M.np,M.np); t0 = 0; r0 = 0;
k  =1;
tic
for tol = [10^2*M.opt.dttol M.opt.dttol]% tol for burn in
    
    for l = 1:M.np   % cycle through populations
        M.SS(l).SS(:,1) = M.SS(l).SS(:,k);
        M.GV(l,:,1) = M.GV(l,:,k);
    end
        
    k  =1;
    r =  Inf*ones(M.opt.N,1);
    dt = M.opt.dt;
    M.opt.t = zeros(M.opt.N,1);
    M.opt.Drawplots =  M.opt.Drawplots +1;
    
    while k < M.opt.N && r(k)>10*tol % integrate through time
        
        % --- Prediction step ---
        
        k=k+1;
        for l = 1:M.nc % Compute recieved conductances
            Gout(l,:) = P.A(l).M*squeeze(M.GV(:,l,k-1)) + P.C(l).M;
        end
        
        for l = 1:M.np   % cycle through populations
            
            P1      = M.P(l);
            SSS2    = [reshape(M.SS(l).SS(:,k-1),[],1,1); squeeze(M.GV(l,:,k-1))'];
            [Qdpdt] = fx_LIFpopME(Gout(:,l),P1);
            M1      = expm(.5*dt*Qdpdt);
            SSS     = M1*SSS2; % Midpoint
            SSS1    = M1*SSS;  % Endpoint
            
            M.SS(l).SS(:,k) = SSS(1:P1.LVV);
            M.SS0(l).SS = SSS1(1:P1.LVV);
            
            M.GV(l,:,k)     = SSS(P1.LVV+1:end);
            M.GV0(l,:)      = SSS1(P1.LVV+1:end);
            
        end
        
        % --- Correction step ---
        
        
        for l = 1:M.nc % Compute recieved conductances
            Gout(l,:) = P.A(l).M*(squeeze(M.GV(:,l,k))) + P.C(l).M; % Midpoint
        end
        
        for l = 1:M.np   % cycle through populationsexp(-6)
            P1      = M.P(l);
            SSS2    = [reshape(M.SS(l).SS(:,k-1),[],1,1); squeeze(M.GV(l,:,k-1))'];
            [Qdpdt] = fx_LIFpopME(Gout(:,l),P1);
            SSS     = expm(dt*Qdpdt)*SSS2;
            
            M.SS(l).SS(:,k) = SSS(1:P1.LVV);
            M.GV(l,:,k)     = SSS(P1.LVV+1:end);
        end
        
        M.opt.t(k) = M.opt.t(k-1) + dt;
        
        D = (M.GV0 - squeeze(M.GV(:,:,k))).^2;
        D = sum(D(:));
        
        for l=1:M.np
            D = D + sum((M.SS0(l).SS - squeeze(M.SS(l).SS(:,k))).^2);
        end
        
        
        
        if D >tol
            dt = dt/2;
            k  = k - 1;
%             disp(['timestep -: ' num2str(dt)])
        elseif (sum(D(:)) < tol/8) %&& dt<M.opt.dt*2
            dt = 2*dt;
            inlineprogress(k,M.opt.N)
%             disp(['timestep +: ' num2str(dt)])
        else
            inlineprogress(k,M.opt.N)
        end
        
        if ~(sum(D(:)) > tol)
            
            if M.opt.t(k) > .005
                [r0, t0] = convergence_check(M.GV,M.opt.t,k,M.opt.T0);
                r(k) = r0;
            end
            % -----------------------
            % --- LFP computation ---
            % -----------------------
            
            if ~isfield(M.opt,'wLFP'); M.opt.wLFP = ones(size(M.P(1).T)); end;
            if ~isfield(M.opt,'LocalElectrode'); M.opt.LocalElectrode = 0; end;
            
            for l = 1:M.np
                P1 = M.P(l);
                FV = P1.FvarV;
                FV(P1.Tr+1:end,:) = FV(repmat(P1.R,P1.LQ,1),:); % heuristic on refractory period
                M.LFP(l,k) = M.opt.wLFP*(Gout(:,l).*(FV'*squeeze(M.SS(l).SS(:,k)*P1.Vres))); % [nS*mV] = [pA] in pico Ampere (per neuron)
                if M.opt.LocalElectrode
                    M.LFP(l,k) = M.LFP(l,k) + squeeze(M.SS(l).SS(P1.Tr+1,k)*P1.Vres)*P1.C*(P1.Vt-P1.Vr)/.001; % Current for repolarization in 1 ms (per neuron)
                end
            end
            
            % --- Current computation ---
            for m= 1:M.np
                P1  =  M.P(m);
                for l = 1:M.nc
                    M.Currents(m,l,k) = Gout(l,m)*(P1.Vg(1)-P1.Vg(2)) + P1.gl*(P1.Vl-P1.Vg(2)); % [nS*mV] = [pA] in pico Ampere (per neuron)
                end
            end
            
            
            if ~isfield(M.opt,'pop'); M.opt.pop = 2; end;
            pop = M.opt.pop;
            
            if M.opt.Drawplots
                
                for l = 1:M.np   % cycle through populationsexp(-6)
                    P1 = M.P(l);
                    subplot(M.np+2,1,l)
                    plot(P1.VV,(M.SS(l).SS(:,k)));axis([P1.VV(1) P1.VV(end) -.001 .2]);
                    drawnow;
                end
                
                
                subplot(M.np+2,1,M.np+1)
                plot(M.opt.t,M.LFP(pop,:)','k');
                hold on
                plot(M.opt.t,squeeze(M.Currents(pop,:,:))','.','MarkerSize',1)
                hold off
                
                subplot(M.np+2,1,M.np+2)
%                 plot(squeeze(M.Currents(pop,1,:)),squeeze(M.Currents(pop,2,:)),'.','MarkerSize',1);
                plot(M.opt.t,r,'.'); hold on;plot(t0,r0,'r.'); hold off;
                xlabel('E (pA)'); %xlim([for l = 1:M.np
                ylabel('I (pA)'); %ylim([-100 100]+EIcurrents(k,2));
                drawnow;
            end
        end
        
        
    end
    
end

toc

Kend = k;
% ---- Plots --------------------------------------------------------------

% for l = 1:M.np   % cycle through populations
%     P1 = M.P(l);
%     subplot(M.np+2,1,l)
%     plot(P1.VV,(M.SS(l).SS(:,k)));axis([P1.VV(1) P1.VV(end) -.001 .2]);
% end
% 
% subplot(M.np+2,1,M.np+1)
% plot(dt:dt:N*dt,LFP(pop,:),'k');
% hold on
% plot(dt:dt:N*dt,EIcurrents)
% hold off
% 
% subplot(M.np+2,1,M.np+2)
% plot(EIcurrents(:,1),EIcurrents(:,2));
% xlabel('E (pA)'); %xlim([-100 100]+EIcurrents(k,1));
% ylabel('I (pA)'); %ylim([-100 100]+EIcurrents(k,2));
% drawnow;
% disp(k*dt)

% ----\Plots\--------------------------------------------------------------

%% dxdp
T0  = M.opt.t(Kend)-t0;
t   = M.opt.t;
NT0 = ceil(2.5*sum(t>t0));

M.GVu        = zeros(M.np,M.nc,NT0);
M.GVu(:,:,1) = M.GV(:,:,Kend);
M.Currentsu  = zeros(M.np,M.nc,NT0);
M.GVu0       = zeros(M.np,M.nc);

for k = 1:M.np
    [U S V]      = svd(M.SS(k).SS(:,1:Kend));
    cs           = cumsum(diag(S));
    no           = sum(cs<(1- M.opt.svdtol)*cs(end));
    M.SS(k).SSu  = [U(:,1:no),      zeros(size(U,1),M.nc);
                    zeros(M.nc,no), eye(M.nc)];
    M.SS(k).SSuT = zeros(no,NT0);
   
    M.SS(k).SSuT(:,1) = U(:,1:no)\M.SS(k).SS(:,Kend);
end



Vp1      = [Vp zeros(size(Vp,1),1)];
dxdp     = zeros([size(M.Currentsu) size(Vp1,2)]);
M.opt.tu = zeros(NT0,size(Vp1,2));
M.fu     = zeros(1,size(Vp1,2));

s          = size(dxdp);
s(3)       = M.opt.Ncoefs;
M.LCrephased = zeros(s);

for j = 1:size(Vp1,2)
    
    disp(['Computing gradient in direction: ' num2str(j)])
    
    P  = spm_ExpP(spm_unvec(spm_vec(P0) +  M.opt.dpsize*Vp1(:,j),P0));
    M  = spm_lifpopsys_LC_adapt(P,M);
    dt = M.opt.dt;
    k  = 1;    
    while k < NT0 % integrate through time
        
        k = k+1;
        % --- Prediction step ---dUp\Qdpdt*Up
        
        %     k=k+1;
        for l = 1:M.nc % Compute recieved conductances
            Gout(l,:) = P.A(l).M*squeeze(M.GVu(:,l,k-1)) + P.C(l).M;
        end
        
        for l = 1:M.np   % cycle through populations
            
            P1      = M.P(l);
            Up      = M.SS(l).SSu;
            SSS2    = [reshape(M.SS(l).SSuT(:,k-1),[],1,1); squeeze(M.GVu(l,:,k-1))'];
            [Qdpdt] = fx_LIFpopME(Gout(:,l),P1);
            M1      = expm(.5*dt*(Up\Qdpdt*Up));
            SSS     = M1*SSS2; % Midpoint
            SSS1    = M1*SSS;  % Endpoint
            
            M.SS(l).SSuT(:,k) = SSS(1:end-M.nc);
            M.GVu(l,:,k)      = SSS(end-M.nc+1:end);
            
            M.SS0(l).SSuT   = SSS1(1:end-M.nc);
            M.GVu0(l,:)     = SSS1(end-M.nc+1:end);
            
        end
        
        % --- Correction step ---
        
        for l = 1:M.nc % Compute recieved conductances
            Gout(l,:) = P.A(l).M*(squeeze(M.GVu(:,l,k))) + P.C(l).M; % Midpoint
        end
        
        for l = 1:M.np   % cycle through populations
            P1      = M.P(l);
            Up      = M.SS(l).SSu;
            SSS2    = [reshape(M.SS(l).SSuT(:,k-1),[],1,1); squeeze(M.GVu(l,:,k-1))'];
            [Qdpdt] = fx_LIFpopME(Gout(:,l),P1);
            SSS     =  expm(dt*(Up\Qdpdt*Up))*SSS2;
            
            M.SS(l).SSuT(:,k) = SSS(1:end-M.nc);
            M.GVu(l,:,k)     = SSS(end-M.nc+1:end);
        end
        
        M.opt.tu(k,j) = M.opt.tu(k-1,j) + dt;
        
        Du = (M.GVu0 - squeeze(M.GVu(:,:,k))).^2;
        Du = sum(Du(:));
        
        for l=1:M.np
            Du = Du + sum((M.SS0(l).SSuT - squeeze(M.SS(l).SSuT(:,k))).^2);
        end
        
        
        
        if Du > M.opt.dttol
            dt = dt/2;
            k  = k - 1;
%             disp(['timestep -: ' num2str(dt)])
        elseif (sum(Du(:)) < M.opt.dttol/8) %&& dt<M.opt.dt*2
            dt = 2*dt;
            inlineprogress(k,NT0)
%             disp(['timestep +: ' num2str(dt)])
        else
            inlineprogress(k,NT0)
        end
        
        if ~(sum(Du(:)) > M.opt.dttol)
            
            for m= 1:M.np
                P1  =  M.P(m);
                for l = 1:M.nc
                    M.Currentsu(m,l,k) = Gout(l,m)*(P1.Vg(1)-P1.Vg(2)) + P1.gl*(P1.Vl-P1.Vg(2)); % [nS*mV] = [pA] in pico Ampere (per neuron)
                end
            end
        end
        
    end
    
    dxdp(:,:,:,j) = M.Currentsu;
    
    [r0, t0] = convergence_check(M.Currentsu,M.opt.tu(:,j),NT0,1.5*T0);
    M.fu(j)  = 1/(M.opt.tu(end,j) - t0);
    M.LCrephased(:,:,:,j) = spm_rephase(squeeze(dxdp(:,:,:,j)),M.opt.tu(:,j),M.opt.Ncoefs,M.fu(j),P.ang);
    
end


%%

x.lc =  shiftdim(M.LCrephased(:,:,:,end),2);
x.f  =  M.fu(end);

dxdp.lc = bsxfun(@minus,x.lc,shiftdim(M.LCrephased(:,:,:,1:end-1),2));
dxdp.f  = M.fu(1:end-1) - x.f;



%%
% 
% figure;
% plot(squeeze(M.GV(1,1,Kend-floor(NT0/2.5):Kend))',squeeze(M.GV(2,1,Kend-floor(NT0/2.5):Kend))','r')
% hold on
% plot(squeeze(M.GVu(1,1,:))',squeeze(M.GVu(2,1,:))')
% 
% %%
% 
% figure;
% plot(squeeze(x(1,1,:))',squeeze(x(1,2,:))','b')
% hold on
% plot(squeeze(dxdp(1,1,:,1))',squeeze(dxdp(1,2,:,1))','.r')
% 
%%

% figure;
% plot(squeeze(M.LCrephased(1,2,:,2))-squeeze(M.LCrephased(1,2,:,1)),'r')


% %%
% 
% for l = 1:M.np
% 
%         M.SS(l).SS(:,1) = M.SS(l).SS(:,end);
%         M.GV(l,:,1)     = M.GV(l,:,end);
% 
% end
% 

