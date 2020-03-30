function [f,J,Q] = spm_fx_bg(x,u,P,M)

% state equations for a neural mass model of the cortex - including E & I populations
% and also the recurrent STN-GPe circuit using Wilson-Cowan type equations

% order               states
% 1 = STN             x(1,1)
% 2 = GPE             x(1,2)
% 3 = CTX 1-E         x(1,3)
% 4 = CTX 2-I         x(1,4)

% delays between structures for subsequent steps
p.dTsg = 6e-3;   % s
p.dTgs = 6e-3;   % s
p.dTgg = 4e-3;   % s
p.dTcs = 6e-3; % s 
p.dTsc = 25e-3;  % s
p.dTcc = 5e-3;   % s

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
    D = [1e-3     p.dTsg  p.dTsc 1e-3      ;...
         p.dTgs   p.dTgg  1e-3   1e-3      ;...
         p.dTcs   1e-3    1e-3   p.dTcc ;...
         1e-3     1e-3    p.dTcc p.dTcc ];
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
p.Te   = 10e-3;  % ms
p.Ti   = 10e-3;  % ms


% BASELINE FIRING RATES
p.Ctx  = 172.18; %  27;     % spm/s
p.Str  = 8.46;   %  2;      % sp/s
p.Ms   = 300;    % sp/s
p.Bs   = 17;     % sp/s
p.Mg   = 400;    % sp/s
p.Bg   = 75;     % sp/s
p.Mi   = 205;    % normal range 200-330 sp/s
p.Bi   = 9.87;   % normal range 0-20 sp/s
p.Me   = 75;     % normal range 50-80 sp/s
p.Be   = 17.85;  % normal range 0-20 sp/s

% for the cortical populations
p.dTcs  = 6e-3;% ms - allowed range 1-10 ms
p.dTcc  = 4e-3;  % ms - allowed range 1-10 ms

p.Wsg  = 19;%4.87;
p.Wgs  = 9;%1.33;
p.Wgg  = 6;%0.53;
p.Wcs  = 9;

% cortical parameters
p.Wsc  = 0.4;
p.Wcc  = 3;

% Y(1) = STN
% Y(2) = GP
% Y(3) = Cortex - Excitatory
% Y(4) = Cortex - Inhibitory

% Add some noise or inputs to fixed values
p.Str     = p.Str + U;
p.Ctx     = p.Ctx + U;

dfdt      =    sparse(n,4);
dfdt(:,1) =   (sigmoid(-p.Wgs* x(:,2) + p.Wcs* x(:,3),p.Ms,p.Bs) - x(:,1))/p.Ts;          % STN
dfdt(:,2) =   (sigmoid( p.Wsg* x(:,1) - p.Wgg* x(:,2) - p.Str,p.Mg,p.Bg) - x(:,2))/p.Tg;  % GP
dfdt(:,3) =   (sigmoid(-p.Wsc* x(:,1) - p.Wcc* x(:,4) + p.Ctx,p.Me,p.Be) - x(:,3))/p.Te;  % CTX 1
dfdt(:,4) =   (sigmoid(p.Wcc*  x(:,3) - p.Wcc* x(:,4),p.Mi,p.Bi) - x(:,4))/p.Ti;          % CTX 2

f       = dfdt'; 

if nargout<2; return,end

% Jacobian for delays (evaluate numerically for simplicity)
%==========================================================================
dfdx = spm_diff('spm_fx_bg',x,u,P,M,1);
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
a     = 1/4;
TOL   = norm(QJ,Inf)*1e-6;
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
