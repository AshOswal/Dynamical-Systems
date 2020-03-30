function [flag ] = spm_DCM_lifpopsys_LC_adapstep_y_par_PRC(P,M,Vp)
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
% Copyright (C) 2014

%% Initialization for Integration -----------------------------------------
% -------------------------------------------------------------------------
try P = P0; end

if ~spm_lifpopsys_LC_check(P,M,Vp); return; end; % check inputs

P0  = P;                % Save original parameters
P   = spm_ExpP(P0);
M   = spm_lifpopsys_LC_prepare(P,M);



% Integrate over time ----------------------------------------------------
% -------------------------------------------------------------------------

Gout       = zeros(M.nc,M.np);
M.LFP      = zeros(M.np,M.opt.N);
% M.Currents = zeros(M.np,M.nc,M.opt.N);
M.Currents = zeros(size(M.opt.popCurrents,1),M.opt.N);

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
if ~isfield(M.opt,'LCtol'); M.opt.LCtol = 10^-3; end;

M.Gv0 = zeros(M.np,M.np); t0 = 0; r0 = 0;
k  =1;

fprintf('\n Burn in:')
for tol = [10^2*M.opt.dttol M.opt.dttol]% tol for burn in
    
    for l = 1:M.np   % cycle through populations
        M.SS(l).SS(:,1) = M.SS(l).SS(:,k);
        M.GV(l,:,1) = M.GV(l,:,k);
    end
        
    k  =1;
    r =  Inf*ones(M.opt.N,1);
    dt = M.opt.dt;
    M.opt.t = zeros(M.opt.N,1);
    
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
            
            if M.opt.t(k) > .005 && ~mod(k,10)
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

            for m= 1:size(M.opt.popCurrents,1)
                P1              = M.P(M.opt.popCurrents(m,1));
                M.Currents(m,k) = - P1.gl*(P1.Vl - M.opt.popCurrents(m,2)); % [nS*mV] = [pA] in pico Ampere (per neuron)
                for l = 1:M.nc
                    M.Currents(m,k) = M.Currents(m,k) - Gout(l,M.opt.popCurrents(m,1))*(P1.Vg(l) - M.opt.popCurrents(m,2)); % [nS*mV] = [pA] in pico Ampere (per neuron)
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
                plot(M.opt.t(1:k),M.LFP(pop,1:k)','k');
                hold on
                plot(M.opt.t,M.Currents','.','MarkerSize',1)
                hold off
                
                subplot(M.np+2,1,M.np+2)
                plot(- M.Currents(2,2:k), M.Currents(1,2:k),'.','MarkerSize',1);
                xlabel('-E (pA)'); %xlim([for l = 1:M.np
                ylabel('I (pA)'); %ylim([-100 100]+EIcurrents(k,2));
                drawnow;
            end
        end
        
        
    end
    
end

Kend = k;


if Kend == M.opt.N || ~checkLC(M.Currents,M.opt.t,Kend,t0,r0,M.opt.LCtol)
    
    disp('Did not reach limit cycle!')
    flag       = 1;
    Vp1        = [Vp zeros(size(Vp,1),1)];
    fu         = zeros(1,size(Vp1,2));
    s          = [size(M.Currents) size(Vp1,2)];
    s(2)       = M.opt.Ncoefs;
    LCrephased = zeros(s);
    FR         = zeros(size(M.opt.popFiringRates,1), size(Vp1,2));
    
%     x.lc =  NaN*permute(LCrephased(:,:,end),[2 1]);
%     x.f  =  NaN*fu(end)*ones(M.opt.Ncoefs,1);
%     
%     dxdp.lc = NaN*permute(bsxfun(@minus,LCrephased(:,:,1:end-1),LCrephased(:,:,end)),[2 1 3]);
%     dxdp.f  = NaN*ones(M.opt.Ncoefs,1)*fu(1:end-1) -  fu(end);
%     dxdp    = cat(2,dxdp.lc,reshape(dxdp.f,[size(dxdp.f,1) 1 size(dxdp.f,2)]))./M.opt.dpsize;

    x.lc =  NaN*permute(LCrephased(:,:,end),[2 1]);
    x.f  =  NaN*fu(end)*ones(M.opt.Ncoefs,1);
    x.fr =  NaN*ones(M.opt.Ncoefs,1)*FR(:,end)';
    
    dxdp.lc = NaN*permute(bsxfun(@minus,LCrephased(:,:,1:end-1),LCrephased(:,:,end)),[2 1 3]);
    dxdp.f  = NaN*ones(M.opt.Ncoefs,1)*fu(1:end-1) -  fu(end);
    dxdp.f  = reshape(dxdp.f,[size(dxdp.f,1) 1 size(dxdp.f,2)]);
    dxdp.fr = NaN*bsxfun(@minus,FR(:,1:end-1),FR(:,end));
    dxdp.fr = repmat(reshape(dxdp.fr,[1 size(dxdp.fr)]),[M.opt.Ncoefs 1 1]);
    dxdp    = cat(2,dxdp.lc,dxdp.f,dxdp.fr)./M.opt.dpsize;
    
    return
else
    flag = 0;
end


%% dxdp
mindt = min(diff(M.opt.t(1:Kend)));
T0  = M.opt.t(Kend)-t0;
t   = 0:mindt:4*T0;
NT0 = length(t);%ceil(2.5*sum(t>t0));

M.GVu        = zeros(M.np,M.nc,NT0);
M.GVu(:,:,1) = M.GV(:,:,Kend);
% M.Currentsu  = zeros(M.np,M.nc,NT0);
M.Currentsu  = zeros(size(M.opt.popCurrents,1),NT0);
M.GVu0       = zeros(M.np,M.nc);

M.opt.svdtol = 0; % FOR PRC - screw efficiency

for k = 1:M.np
    [U, S, ~]      = svd(M.SS(k).SS(:,1:Kend)); % U = eye(size(U)); S = U;
    cs           = cumsum(diag(S));
    no           = sum(cs<=(1-M.opt.svdtol)*cs(end));
    M.SS(k).SSu  = [U(:,1:no),      zeros(size(U,1),M.nc);
                    zeros(M.nc,no), eye(M.nc)];
    M.SS(k).SSuT = zeros(no,NT0);
    M.SS(k).SSuT(:,1) = U(:,1:no)\M.SS(k).SS(:,Kend);
end



Vp1      = [Vp zeros(size(Vp,1),1)];
M.opt.tu = zeros(NT0,10);
fu       = zeros(1,10);

s          = [size(M.Currentsu) 10];
s(2)       = M.opt.Ncoefs;
LCrephased = zeros(s);
FR         = zeros(size(M.opt.popFiringRates,1), 10);

fprintf('\n Computing Gradients:')
ndp = 10;

mkdir PRCi

%%

jmax = 50;

for l = 1:10
    
for j = 1:jmax
    disp([l, j])
    Pj  = spm_ExpP(spm_unvec(spm_vec(P0),P0));
    Mj  = spm_lifpopsys_LC_adapt(Pj,M);
    
    k  = 1;
    
    U = 1:NT0;
    U = exp((l-5)/2)*exp(-(1/.002)*(U*mindt-j*T0/jmax)).*(U*mindt > j*T0/jmax);
    
    while k < NT0 % integrate through time
        dt =mindt;
        k = k+1;
        % --- Prediction step ---dUp\Qdpdt*Up
        
        %     k=k+1;
        Goutj = zeros(Mj.nc,Mj.np);
        for l1 = 1:Mj.nc % Compute recieved conductances
            Goutj(l1,:) = Pj.A(l1).M*squeeze(Mj.GVu(:,l1,k-1)) + Pj.C(l1).M;
        end
        
%         Goutj(1,2) = Goutj(1,2) + U(k);
        Goutj(1,1) = Goutj(1,1) + U(k);
        
        for l1 = 1:Mj.np   % cycle through populations
            
            P1      = Mj.P(l1);
            Up      = Mj.SS(l1).SSu;
            SSS2    = [reshape(Mj.SS(l1).SSuT(:,k-1),[],1,1); squeeze(Mj.GVu(l1,:,k-1))'];
            [Qdpdt] = fx_LIFpopME(Goutj(:,l1),P1);
            M1      = expm(.5*dt*(Up\Qdpdt*Up));
            SSS     = M1*SSS2; % Midpoint
            SSS1    = M1*SSS;  % Endpoint
            
            Mj.SS(l1).SSuT(:,k) = SSS(1:end-Mj.nc);
            Mj.GVu(l1,:,k)      = SSS(end-Mj.nc+1:end);
            
            Mj.SS0(l1).SSuT   = SSS1(1:end-Mj.nc);
            Mj.GVu0(l1,:)     = SSS1(end-Mj.nc+1:end);
            
        end
        
        % --- Correction step ---
        
        for l1 = 1:Mj.nc % Compute recieved conductances
            Goutj(l1,:) = Pj.A(l1).M*(squeeze(Mj.GVu(:,l1,k))) + Pj.C(l1).M; % Midpoint
        end
        
%         Goutj(1,2) = Goutj(1,2) + U(k);
        Goutj(1,1) = Goutj(1,1) + U(k);
        
        for l1 = 1:Mj.np   % cycle through populations
            P1      = Mj.P(l1);
            Up      = Mj.SS(l1).SSu;
            SSS2    = [reshape(Mj.SS(l1).SSuT(:,k-1),[],1,1); squeeze(Mj.GVu(l1,:,k-1))'];
            [Qdpdt] = fx_LIFpopME(Goutj(:,l1),P1);
            SSS     =  expm(dt*(Up\Qdpdt*Up))*SSS2;
            
            Mj.SS(l1).SSuT(:,k) = SSS(1:end-Mj.nc);
            Mj.GVu(l1,:,k)     = SSS(end-Mj.nc+1:end);
        end
        
%         Mj.opt.tu(k,j) = Mj.opt.tu(k-1,j) + dt;
        
%         Du = (Mj.GVu0 - squeeze(Mj.GVu(:,:,k))).^2;
%         Du = sum(Du(:));
%         
%         for l1=1:Mj.np
%             Du = Du + sum((Mj.SS0(l1).SSuT - squeeze(Mj.SS(l1).SSuT(:,k))).^2);
%         end
        
%         if Du > M.opt.dttol && firstiter
%             dt = dt/2;
%             k  = k - 1;
%         elseif (sum(Du(:)) < M.opt.dttol/8) && firstiter 
%             dt = 2*dt;
%         elseif ~firstiter && k<NT0
%             dt =  M.opt.tu(k+1,end) -  M.opt.tu(k,end);
%         end
        
        % --- Compute currents ---
        % ------------------------
        if 1 %~(sum(Du(:)) > Mj.opt.dttol && firstiter)
            for m= 1:size(Mj.opt.popCurrents,1)
                P1              = Mj.P(Mj.opt.popCurrents(m,1));
                Mj.Currentsu(m,k) = - P1.gl*(P1.Vl - Mj.opt.popCurrents(m,2)); % [nS*mV] = [pA] in pico Ampere (per neuron)
                for l1 = 1:Mj.nc
                    Mj.Currentsu(m,k) = Mj.Currentsu(m,k) - Goutj(l1,Mj.opt.popCurrents(m,1))*(P1.Vg(l1) - Mj.opt.popCurrents(m,2)); % [nS*mV] = [pA] in pico Ampere (per neuron)
                end
            end
        end
        
    end
    
%     [~, t1] = convergence_checkC(Mj.Currentsu,Mj.opt.tu(:,j),NT0,1.5*T0);
%     fu(j)  = 1/(Mj.opt.tu(end,j) - t1);
%     LCrephased(:,:,j) = spm_rephase(Mj.Currentsu,Mj.opt.tu(:,j),Mj.opt.Ncoefs,fu(j),Pj.ang);
%     FR(:,j) = mean(spm_rephase(squeeze(Mj.GVu(M.opt.popFiringRates,1,:)),Mj.opt.tu(:,j),Mj.opt.Ncoefs,fu(j),Pj.ang),2);
       
    fprintf('.')
    
%     s0 = 0;
%     s1 = 0;
%     for l1 = 1:Mj.np
%         s0 = s0 + sum(Mj.SS(l1).SSu(1:end - Mj.nc,1:end - Mj.nc)*Mj.SS(l1).SSuT(:,end));
%         s1 = s1 + sum(Mj.SS(l1).SSu(1:end - Mj.nc,1:end - Mj.nc)*Mj.SS(l1).SSuT(:,1));
%     end
%     if abs(s0-s1) > .01*s1
%         disp(['Warning: M.opt.svdtol too high! ratio:' num2str(s0/s1)])
%     end
    save(['PRCi/Mj' num2str(l) '_' num2str(j)], 'Mj')
    
end

end

%

% disp(['True freq: ' num2str(1/T0) ' SVD freq: ' num2str(x.f(1))])
% % 
% figure;plot(squeeze(dxdp(:,1,:)))
% figure;plot((squeeze(M.LCrephased(1,:,:))))




