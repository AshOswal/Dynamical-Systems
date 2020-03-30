function y = spm_fx_if(x)

% state equation for IF unit
% A Oswal

% x(1) = input - poisson spikes  
% x(2) = transmembrane potential
% x(3) = time since last spike

% fixed parameters
%--------------------------------------------------------------------------
C     =   0.375;                % Capacitance {nF}
Vl    =     -73;                % Resting potential {mV}
gl    =      25;                % passive conductance {nS}
Vd    =     -63;                % spiking voltage {mV}
Rm    =      100;                 % resistance 
gc    =   (Vd - Vl)/Rm;
input =   sum(x(1))*gc*Rm;

% dV/dt {mV/s)
%--------------------------------------------------------------------------
v     = x(2);
T     = x(3);
dVdt  = (1/C) * (gl*(Vl   - v)) + input;         
s     = exp(-T^2/(2*(1e-3)^2)); % a gaussian blur cf. refractory
dVdt  = dVdt + 1e4*(-90 - v)*s;


% pst
%--------------------------------------------------------------------------
dtdt  = 1 - T*1e4*(v > Vd);

% dx/dt
%--------------------------------------------------------------------------
y     = [dVdt; dtdt];