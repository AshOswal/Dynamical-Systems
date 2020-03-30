function M   = spm_lifpopsys_LC_prepare(P,M)



% if isfield(M,'P'); return; end
if ~isfield(M.opt,'r'); M.opt.r = 1; end

M.P     = P.P;
M.np    = length(P.P);      % number o populations
M.nc    = length(P.P(1).T); % number of channels (ion)

%% initialization of populations ---------------------------------

for k = 1:M.np
    
%     M.P(k).Vt  = -42;            % [mV] Spike threshold voltage
%     M.P(k).Vr  = -70;            % [mV] Reset voltage
    M.P(k).d   = 0/1000;         % [s]  Afferent axonal delay
    M.P(k).G   = 0*P.P(k).T + 1;         % [1]  Synaptic eficacy scaling (not in use)
    M.P(k).DL  = 0*P.P(k).T;       % Unused parameters to model multiplicative membrane noise
    M.P(k).DGL = 0*P.P(k).T;       % Unused parameters to model multiplicative membrane noise
    
    M.P(k).Vmin = -80;           % [mV] minimum of the voltage space
    M.P(k).Tr   = M.opt.r*80;            % number of bins for Voltage space discretization
    M.P(k).LQ   = M.opt.r*20;            % number of bins for Queue space discretization
    
    % -------- Do not Change --------
    M.P(k).LVV  = M.P(k).Tr+M.P(k).LQ;
    M.P(k).VV   = linspace(M.P(k).Vmin,M.P(k).Vt,M.P(k).Tr); % voltage discretisation
    M.P(k).Vres = M.P(k).VV(2)-M.P(k).VV(1);                 % interval
    M.P(k).VV   = [M.P(k).VV (M.P(k).VV(end)+M.P(k).Vres):M.P(k).Vres:(M.P(k).VV(end)+M.P(k).LQ*M.P(k).Vres)];
    Fsteady = zeros(M.P(k).LVV,1);
    FvarV   = zeros(M.P(k).LVV,length(M.P(k).G));
    [C R]   = min(abs(M.P(k).VV-M.P(k).Vr)); % C - min value R - index
    M.P(k).R     = R; % reset bin - a small discretization error is made here
    for l = 1:M.P(k).Tr % no. bins voltage discretisation
        Fsteady(l) = 1/M.P(k).C*(-M.P(k).gl*(M.P(k).VV(l)-M.P(k).Vl)) ;
        for m = 1:length(M.P(k).G)
            FvarV(l,m) = 1/M.P(k).C*(M.P(k).VV(l)-M.P(k).Vg(m)); % this is as above but multiplied by the time dependent conductances later
        end
    end
    Fsteady(M.P(k).Tr+1:end) = M.P(k).LQ/M.P(k).Rt; % during the queue period
    M.P(k).delay = min(M.P(k).Tr + 1 + floor(M.P(k).LQ*M.P(k).d/M.P(k).Rt),M.P(k).LVV-1);
    M.P(k).FvarV   = FvarV;
    M.P(k).Fsteady = Fsteady;

end

%% Initialization for Integration -----------------------------------------

if ~isfield(M,'opt'); M.opt = []; end;
if ~isfield(M.opt,'N'); M.opt.N = 1000; end;

if isfield(M,'SS')
    for k = 1:M.np
        M.SS(k).SS(:,1) = M.SS(k).SS(:,end);
    end
    M.GV(:,:,1) = M.GV(:,:,end);
else
    M.GV = zeros(M.np,M.nc);
    for k = 1:M.np
        M.SS(k).SS  = null(fx_LIFpopME(0*M.P(k).T,M.P(k))); % initializes populations with stationary probability density functions (pdfs) - dpdt = 0.
        M.SS(k).SS  = M.SS(k).SS/sum(M.SS(k).SS(1:end-M.nc))/M.P(k).Vres; % normalization to probability density per milivolt
        M.GV(k,:)   = M.SS(k).SS(end-M.nc+1:end)*0;
        M.SS(k).SS  = repmat(M.SS(k).SS(1:end-M.nc),1,M.opt.N);
    end
    M.GV = repmat(M.GV,[1,1,M.opt.N]);
end

%% Check options

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



end

