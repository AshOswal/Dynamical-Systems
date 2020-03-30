function [y] = LIF_StSolLog(Go,P,flag)

% Go = [0.25 0];

A0 = 1/P.C*(P.gl)*[-1 P.Vl]';

Avar = zeros(length(P.G),2);
for m = 1:length(P.G)
    Avar(m,:) = 1/P.C*P.G(m)*[-1 P.Vg(m)];
end

AF = A0 + Avar'*Go;
AF = AF.*[1/P.D; 1/P.D]; %review calculations

Gx  = @(x) -(sqrt(-AF(1))*(x + AF(2)/AF(1))).^2;
% dGx = @(x) -2*(-AF(1)*(x + AF(2)/AF(1))).*exp(-(sqrt(-AF(1))*(x + AF(2)/AF(1))).^2);

Fx  = @(x) DawsonIntLog((sqrt(-AF(1))*(x + AF(2)/AF(1))),flag);  
% dFx = @(x) sqrt(-AF(1))*(1-2*(sqrt(-AF(1))*(x + AF(2)/AF(1))).*DawsonInt((sqrt(-AF(1))*(x + AF(2)/AF(1))),flag));

N1 = 1;
S2 = 1/(-exp(Fx(P.Vr)-Gx(P.Vr)) + exp(Fx(P.Vt)-Gx(P.Vt)));
N2 = Fx(P.Vt)/Gx(P.Vt)*S2;
S3 = -(N2*dGx(P.Vt)-S2*dFx(P.Vt))*P.D/2; % Check the 1/2 ratio!
F3 = (P.Vmax-P.Vt)/P.Rt;


G = Gx(P.VV(1:P.Tr));
F = Fx(P.VV(1:P.Tr));
y = [N1*G(1:P.R) N2*G(P.R+1:end)-S2*F(P.R+1:end) S3/F3*ones(1,P.LVV-P.Tr)];
Sy = sum(y)/P.Vres;
y  = y/Sy;

% plot(y)
% plot([F' G'])
end