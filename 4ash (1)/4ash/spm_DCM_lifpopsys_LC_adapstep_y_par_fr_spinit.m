function [ dxdp, x ] = spm_DCM_lifpopsys_LC_adapstep_y_par_fr_spinit(P,M,Vp)
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

if ~isfield(M.opt,'dttol'); M.opt.dttol = 10^-3; end;
if ~isfield(M.opt,'dt'); M.opt.dt = .001; end;
if ~isfield(M.opt,'T0max'); M.opt.T0max = .12; end; 
if ~isfield(M.opt,'svdtol'); M.opt.svdtol = 10^-6; end; 
if ~isfield(M.opt,'dpsize'); M.opt.dpsize = exp(-3); end;
if ~isfield(M.opt,'Ncoefs'); M.opt.Ncoefs = 1000; end;
if ~isfield(M.opt,'LCtol'); M.opt.LCtol = 1; end;

M.Gv0 = zeros(M.np,M.np); 

M.LFP      = zeros(M.np,M.opt.N);
% M.Currents = zeros(M.np,M.nc,M.opt.N);
M.Currents = zeros(size(M.opt.popCurrents,1),M.opt.N);

if ~isfield(M.opt,'Drawplots')
    M.opt.Drawplots = 0; 
elseif M.opt.Drawplots
    figure(M.opt.Drawplots)
end


%% Burn in


M = spm_DCM_lifpopsys_LC_int_spinit(P,M);



%% dxdp
Kend = M.opt.Kend;
T0  = 1/M.opt.U.f0;
ti  = sum(M.opt.t(1:Kend)<M.opt.t(Kend)-T0)+1;
mindt = min(diff(M.opt.t(ti:Kend))); 
if mindt < 10^-5;  fprintf('bad mindt: %f', mindt); mindt = 10^-5; end
t   = 0:mindt:4*T0;
NT0 = length(t);%ceil(2.5*sum(t>t0));

M.GVu        = zeros(M.np,M.nc,NT0);
M.GVu(:,:,1) = M.GV(:,:,Kend);
M.FiringRu   = zeros(M.np,NT0);
% M.Currentsu  = zeros(M.np,M.nc,NT0);
M.Currentsu  = zeros(size(M.opt.popCurrents,1),NT0);
M.GVu0       = zeros(M.np,M.nc);
sno = 0; % display the total number of components 

for k = 1:M.np
    [U, S, ~]      = svd(M.SS(k).SS(:,ti:Kend)); % U = eye(size(U)); S = U;
    cs           = cumsum(diag(S));
    no           = sum(cs<=(1-M.opt.svdtol)*cs(end)) +1;
    if no<4, no=4; end
    sno          = sno + no; 
    M.SS(k).SSu  = [U(:,1:no),      zeros(size(U,1),M.nc);
                    zeros(M.nc,no), eye(M.nc)];
    M.SS(k).SSuT = zeros(no,NT0);
    M.SS(k).SSuT(:,1) = U(:,1:no)\M.SS(k).SS(:,Kend);
end



Vp1      = [Vp zeros(size(Vp,1),1)];
M.opt.tu = zeros(NT0,size(Vp1,2));
fu       = zeros(1,size(Vp1,2));

s          = [size(M.Currentsu) size(Vp1,2)];
s(2)       = M.opt.Ncoefs;
LCrephased = zeros(s);
FR         = zeros([size(M.opt.popFiringRates,1), M.opt.Ncoefs, size(Vp1,2)]);
% FR         = zeros(size(M.opt.popFiringRates,1), size(Vp1,2));

fprintf('\n Computing Gradients (%i principal dimensions): ', sno)
ndp = size(Vp1,2);

parfor j = 1:ndp
    
    Pj  = spm_ExpP(spm_unvec(spm_vec(P0) +  M.opt.dpsize*Vp1(:,j),P0));
    Mj  = spm_lifpopsys_LC_adapt(Pj,M);
    
    k  = 1;    
    while k < NT0 % integrate through time
        dt =mindt;
        k = k+1;
                
        Mj.opt.tu(k,j) = Mj.opt.tu(k-1,j) + dt;
        
        % --- Prediction step ---dUp\Qdpdt*Up
        
        %     k=k+1;
        Goutj = zeros(Mj.nc,Mj.np);
        for l1 = 1:Mj.nc % Compute recieved conductances
            Goutj(l1,:) = Pj.C(l1).M; %Pj.A(l1).M*squeeze(Mj.GVu(:,l1,k-1)) + 
            if ~all(isfinite(Goutj(l1,:))), error('improper conductance values'); end
        end
        U  = Mj.opt.U.u(Mj.opt.tu(k,j)); % external inputs
        Goutj = Goutj + U;

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
            Goutj(l1,:) = Pj.C(l1).M; % Pj.A(l1).M*(squeeze(Mj.GVu(:,l1,k))) +; Midpoint
        end
        U  = Mj.opt.U.u(Mj.opt.tu(k,j) + .5*dt);
%         U  = M.opt.U.u(M.opt.t(k) + .5*dt); % external inputs
        Goutj = Goutj + U;
        
        for l1 = 1:Mj.np   % cycle through populations
            P1      = Mj.P(l1);
            Up      = Mj.SS(l1).SSu;
            SSS2    = [reshape(Mj.SS(l1).SSuT(:,k-1),[],1,1); squeeze(Mj.GVu(l1,:,k-1))'];
            [Qdpdt] = fx_LIFpopME(Goutj(:,l1),P1);
            SSS     =  expm(dt*(Up\Qdpdt*Up))*SSS2;
            
            Mj.SS(l1).SSuT(:,k) = SSS(1:end-Mj.nc);
            Mj.GVu(l1,:,k)      = SSS(end-Mj.nc+1:end);
            
            
            Mj.FiringRu(l1,k)  = Up(P1.Tr+1,:)*SSS.*P1.Fsteady(P1.Tr+1).*P1.Vres;
        end

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
            % - Firing rates
%             for m= 1:size(M.opt.popFiringRates,1)
%                 P1              = Mj.P(Mj.opt.popFiringRates(m,1));
%                 Mj.FiringRu(m,k) = - P1.gl*(P1.Vl - Mj.opt.popCurrents(m,2));
%             end
%             
            
        end
        
    end
    
%     [~, t1] = convergence_checkCSVD(Mj.Currentsu,Mj.opt.tu(:,j),NT0,1.5*T0);
    fu(j)   = (Mj.opt.U.f0);%1/(Mj.opt.tu(end,j) - t1);
    
    %rephase currents and firing rates
    Original = [bsxfun(@minus,Mj.Currentsu,mean(Mj.Currentsu,2));
                squeeze(Mj.FiringRu(M.opt.popFiringRates,:,:))];
    Rephased = spm_rephase(Original,Mj.opt.tu(:,j),Mj.opt.Ncoefs,fu(j));
    LCrephased(:,:,j) = Rephased(1:size(Mj.Currentsu,1),:);
    LCrephased(:,:,j) = bsxfun(@minus,LCrephased(:,:,j),mean(LCrephased(:,:,j),2));
    FR(:,:,j)         = Rephased(size(Mj.Currentsu,1)+1:end,:);
    
%     LCrephased(:,:,j) = spm_rephase(Mj.Currentsu,Mj.opt.tu(:,j),Mj.opt.Ncoefs,fu(j),Pj.ang);
%     FR(:,:,j)         = spm_rephase(squeeze(Mj.FiringRu(M.opt.popFiringRates,:,:)),Mj.opt.tu(:,j),Mj.opt.Ncoefs,fu(j),Pj.ang);
%     FR(:,:,j) = (spm_rephase2(squeeze(Mj.GVu(M.opt.popFiringRates,1,:)),Mj.opt.tu(:,j),Mj.opt.Ncoefs,fu(j),Pj.ang));
%     FR(:,j) = mean(spm_rephase2(squeeze(Mj.GVu(M.opt.popFiringRates,1,:)),Mj.opt.tu(:,j),Mj.opt.Ncoefs,fu(j),Pj.ang),2);
       
    fprintf('.')
    
    s0 = 0;
    s1 = 0;
    for l1 = 1:Mj.np
        s0 = s0 + sum(Mj.SS(l1).SSu(1:end - Mj.nc,1:end - Mj.nc)*Mj.SS(l1).SSuT(:,end));
        s1 = s1 + sum(Mj.SS(l1).SSu(1:end - Mj.nc,1:end - Mj.nc)*Mj.SS(l1).SSuT(:,1));
    end
    if abs(s0-s1) > .01*s1
        fprintf('Warning: M.opt.svdtol too high! ratio:%f \n',(s0/s1))
        if abs(s0-s1) > .05*s1
            error('SVD tolerance broken')
        end
    end
    
    
end


%

% x.lc =  permute(LCrephased(:,:,end),[2 1]);
% x.f  =  fu(end)*ones(M.opt.Ncoefs,1);
% x.fr =  ones(M.opt.Ncoefs,1)*FR(:,end)';
x.fr =  permute(FR(:,:,end),[2 1]);
% x.c  =  zeros(M.opt.Ncoefs,1);
% 
% dxdp.lc = permute(bsxfun(@minus,LCrephased(:,:,1:end-1),LCrephased(:,:,end)),[2 1 3]);
% dxdp.f  = ones(M.opt.Ncoefs,1)*fu(1:end-1) -  fu(end);
% dxdp.f  = reshape(dxdp.f,[size(dxdp.f,1) 1 size(dxdp.f,2)]);
% dxdp.fr = bsxfun(@minus,FR(:,1:end-1),FR(:,end));
% dxdp.fr = repmat(reshape(dxdp.fr,[1 size(dxdp.fr)]),[M.opt.Ncoefs 1 1]);
dxdp.fr = permute(bsxfun(@minus,FR(:,:,1:end-1),FR(:,:,end)),[2 1 3])./M.opt.dpsize;

% dxdp.c  = 0*dxdp.f;
% dxdp    = cat(2,dxdp.lc,dxdp.f,dxdp.fr,dxdp.c)./M.opt.dpsize;




% disp(['True freq: ' num2str(1/T0) ' SVD freq: ' num2str(x.f(1))])
% % 
% figure;plot(squeeze(dxdp(:,1,:)))
% figure;plot((squeeze(M.LCrephased(1,:,:))))


end

