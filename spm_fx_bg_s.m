function [f,J,Q] = spm_fx_bg_s(x,u,P,M)

% state equations for a neural mass model of the cortex - including E & I populations
% and also the recurrent STN-GPe circuit using Wilson-Cowan type equations

% order               states
% 1 = STN             x(1,1)
% 2 = GPE             x(1,2)

% delays between structures for subsequent steps
p.dTsg = 6e-3;     % s
p.dTgs = 6e-3;     % s
p.dTgg = 4e-3;     % s
p.dTcs = 5.5e-3;   % s 

% A Oswal

% the parameter descriptions of the original model as per Nevado Holgado 2010
%--------------------------------------------------------------------------

% get dimensions and configure state variables
%--------------------------------------------------------------------------
x  = spm_unvec(x,M.x);               % neuronal states
n  = size(x,1);                      % number of sources


% input connections
%--------------------------------------------------------------------------
if isfield(P,'C')
    C    = P.C;
else
    C    = 1;
end

% input
%==========================================================================
if isfield(M,'u')
    % endogenous input
    %----------------------------------------------------------------------
    U = u(:);
else
    % exogenous input
    %----------------------------------------------------------------------
    U = C*u(:);
end

% delays 
%--------------------------------------------------------------------------
if isfield(M,'D')
    D = M.D;
else
    % note additional self inhibitory feedback to CTX I
    D = [ 0      p.dTsg ;...
         p.dTgs        p.dTgg];...
end
    
% delay approximation
%--------------------------------------------------------------------------
if isfield(M,'Nit')
    N = M.Nit;
else
    N = 2^8;
end


% TIME CONSTANT FOR ACTIVATION
p.Ts   = 6e-3;   % ms
p.Tg   = 14e-3;  % ms


% BASELINE FIRING RATES
p.Ctx  = 27;     % spm/s
p.Str  = 2;      % sp/s
p.Ms   = 300;    % sp/s
p.Bs   = 17;     % sp/s
p.Mg   = 400;    % sp/s
p.Bg   = 75;     % sp/s

p.K = 1;
% disease
p.Wsg  = 19 + p.K*(20-19);
p.Wgs  = 1.12 + p.K*(10.7-1.12);
p.Wgg  = 6.6 + p.K*(12.3 - 6.6);
p.Wcs  = 2.42 + p.K*(9.2 - 2.42);
p.Wxg  = 15.1 + p.K*(139.4 - 15.1);

% Y(1) = STN
% Y(2) = GP

% Add some noise or inputs to fixed values
p.Str     = p.Str + U;
p.Ctx     = p.Ctx + U;

% first with the delays but without the sigmoid - to 1st order
% df_ND_dt      =    sparse(n,2);
% df_ND_dt(:,1) =    (-p.Wgs* x(:,2) + p.Wcs*p.Ctx - x(:,1))/p.Ts;                 % STN
% df_ND_dt(:,2) =    (p.Wsg* x(:,1) - p.Wgg* x(:,2) - p.Wxg*p.Str - x(:,2))/p.Tg;  % GP
% 
% % now the full delayed model
% dfdt          =    sparse(n,2);
% dfdt(:,1)     =    (sigmoid(-p.Wgs* (x(:,2) - p.dTgs*df_ND_dt(:,2)) + p.Wcs*p.Ctx,p.Ms,p.Bs) - x(:,1))/p.Ts;          % STN
% dfdt(:,2)     =    (sigmoid( p.Wsg* (x(:,1) - p.dTsg*df_ND_dt(:,1)) - p.Wgg* (x(:,2) - p.dTgg*df_ND_dt(:,2)) - p.Wxg*p.Str,p.Mg,p.Bg) - x(:,2))/p.Tg;  % GP

dfdt          =    sparse(n,2);
dfdt(:,1)     =    (sigmoid(-p.Wgs* x(:,2) + p.Wcs*p.Ctx,p.Ms,p.Bs) - x(:,1))/p.Ts;          % STN
dfdt(:,2)     =    (sigmoid( p.Wsg* x(:,1) - p.Wgg* x(:,2) - p.Wxg*p.Str,p.Mg,p.Bg) - x(:,2))/p.Tg;  % GP

f       = dfdt'; 

if nargout<2; return,end

% Jacobian for delays (evaluate numerically for simplicity)
%==========================================================================
dfdx = spm_diff('spm_fx_bg_s',x,u,P,M,1);
J = dfdx;


if nargout<3; return,end

% delay operator
%==========================================================================

% delay operator: first-order approximation if N = 0
% Implement: dx(t)/dt = f(x(t - d)) = inv(1 + D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------
if ~N
    Q  = pinv(eye(length(J)) + D.*J);
    return
end

% delay operator: estimated using a Robbins-Monro algorithm
%--------------------------------------------------------------------------
D     = -D;
QJ    = (eye(length(J)) - D.*J)\J;
a     = 1/16;
TOL   = norm(QJ,Inf)*1e-18;
dn    = 1;
Dn    = cell(N,1);
for n = 1:N
    dn    = dn.*D;
    Dn{n} = dn.*J;
end


% initialise (solve matrix polynomial for Q = QJ)
% Q = sum( ((D.^n).*J)*(Q^n)/factorial(n) ): n = 0,1,...
%==========================================================================
for i = 1:N
    
    QJn   = 1;
    Q     = J;
    for n = 1:N
        
        % n-th order Taylor term
        %------------------------------------------------------------------
        QJn           = QJn*QJ/n;
        dQ            = Dn{n}*QJn;
        dQ(isnan(dQ)) = 0;
        Q             = Q + dQ;
        
        % break if convergence
        %------------------------------------------------------------------
        if norm(dQ,Inf) < TOL, break, end
        
    end

    % Robbins-Monro update and break if convergence
    %----------------------------------------------------------------------
    QJ = QJ*(1 - a) + Q*a;
    
    if norm((QJ - Q),Inf) < TOL, break, end
    
end
Q      = Q*spm_inv(J,exp(-16));
%}

function y = sigmoid (x,max,min)
y = max./(1+(exp(-4*x*inv(max))*(max-min)/min));
end

end



