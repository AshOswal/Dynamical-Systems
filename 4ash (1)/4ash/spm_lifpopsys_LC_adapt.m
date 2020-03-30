function M   = spm_lifpopsys_LC_adapt(P,M)

M.P     = P.P;
M.np    = length(P.P);      % number o populations
M.nc    = length(P.P(1).T); % number of channels

%% update populations ---------------------------------

for k = 1:M.np
    
%     M.P(k).Vt  = -42;            % [mV] Spike threshold voltage
%     M.P(k).Vr  = -70;            % [mV] Reset voltage
    M.P(k).d   = 0/1000;         % [s]  Afferent axonal delay
    M.P(k).G   = 0*P.P(k).T + 1;        % [1]  Synaptic eficacy scaling (not in use)
    M.P(k).DL  = 0*P.P(k).T;       % Unused parameters to model multiplicative membrane noise
    M.P(k).DGL = 0*P.P(k).T;       % Unused parameters to model multiplicative membrane noise
    
    M.P(k).Vmin = -80;           % [mV] minimum of the voltage space
    M.P(k).Tr   = M.opt.r*80;        % number of bins for Voltage space discretization
    M.P(k).LQ   = M.opt.r*20;         % number of bins for Queue space discretization
    
    % -------- Do not Change --------
    M.P(k).LVV  = M.P(k).Tr+M.P(k).LQ;
    M.P(k).VV   = linspace(M.P(k).Vmin,M.P(k).Vt,M.P(k).Tr);
    M.P(k).Vres = M.P(k).VV(2)-M.P(k).VV(1);
    M.P(k).VV   = [M.P(k).VV (M.P(k).VV(end)+M.P(k).Vres):M.P(k).Vres:(M.P(k).VV(end)+M.P(k).LQ*M.P(k).Vres)];
    Fsteady = zeros(M.P(k).LVV,1);
    FvarV   = zeros(M.P(k).LVV,length(M.P(k).G));
    [C R]   = min(abs(M.P(k).VV-M.P(k).Vr));
    M.P(k).R     = R; % reset bin - a small discretization error is made here
    for l = 1:M.P(k).Tr
        Fsteady(l) = 1/M.P(k).C*(-M.P(k).gl*(M.P(k).VV(l)-M.P(k).Vl)) ;
        for m = 1:length(M.P(k).G)
            FvarV(l,m) = 1/M.P(k).C*(M.P(k).VV(l)-M.P(k).Vg(m));
        end
    end
    Fsteady(M.P(k).Tr+1:end) = M.P(k).LQ/M.P(k).Rt;
    M.P(k).delay = min(M.P(k).Tr + 1 + floor(M.P(k).LQ*M.P(k).d/M.P(k).Rt),M.P(k).LVV-1);
    M.P(k).FvarV   = FvarV;
    M.P(k).Fsteady = Fsteady;

end

end

