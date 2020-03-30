function [y] = LIF_StSol(Go,P,flag)

% Go = [.25 .0]';

A0 = 1/P.C*(P.gl)*[-1 P.Vl]';

Avar = zeros(length(P.G),2);
for m = 1:length(P.G)
    Avar(m,:) = 1/P.C*P.G(m)*[-1 P.Vg(m)];
end

AF = A0 + Avar'*Go;
AF = AF.*[1/P.D; 1/P.D]; %review calculations -> This is correct becaus P.D already accounts for the 1/2 fatctor

Gx  = @(x) exp(-(sqrt(-AF(1))*(x + AF(2)/AF(1))).^2);
GxL = @(x) (-(sqrt(-AF(1))*(x + AF(2)/AF(1))).^2);
dGx = @(x) -2*(-AF(1)*(x + AF(2)/AF(1))).*exp(-(sqrt(-AF(1))*(x + AF(2)/AF(1))).^2);

Fx  = @(x) DawsonInt((sqrt(-AF(1))*(x + AF(2)/AF(1))),flag);  %change the flag if too slow
dFx = @(x) sqrt(-AF(1))*(1-2*(sqrt(-AF(1))*(x + AF(2)/AF(1))).*DawsonInt((sqrt(-AF(1))*(x + AF(2)/AF(1))),flag)); %change the flag if too slow


% Care is needed against numerical under/overflow

if Fx(P.Vt)/Gx(P.Vt)>.5 %low firing rate
    N2 = 1;
    S2 = Gx(P.Vt)/Fx(P.Vt);
    N1 = (1 - Fx(P.Vr)*Gx(P.Vt)/Gx(P.Vr));
    S3 = -(N2*dGx(P.Vt)-S2*dFx(P.Vt))*P.D/2; % Check the 1/2 ratio!
    F3 = P.LQ/P.Rt;
    
    G  = Gx(P.VV(1:P.Tr));
    F  = Fx(P.VV(1:P.Tr));
    y  = [N1*G(1:P.R) N2*G(P.R+1:end)-S2*F(P.R+1:end) S3/F3*ones(1,P.LVV-P.Tr)];
    Sy = sum(y)/P.Vres;
    y  = y/Sy; y2=y;
else % high firing rate
    S2  = 1;
    N2  = Fx(P.Vt)/Gx(P.Vt);
    N1a = log(-Fx(P.Vr))-GxL(P.Vr);
    N1b = Fx(P.Vt)/Gx(P.Vt);
    S3  = -(N2*dGx(P.Vt)-S2*dFx(P.Vt))*P.D/2; % Check the 1/2 ratio!
    F3  = (P.LQ)/P.Rt;
       
    G  = Gx(P.VV(1:P.Tr));
    F  = Fx(P.VV(1:P.Tr));
    GL = GxL(P.VV(1:P.R));
    y  = [exp(N1a+GL(1:P.R))+N1b*G(1:P.R) N2*G(P.R+1:end)-S2*F(P.R+1:end) S3/F3*ones(1,P.LVV-P.Tr)];
    Sy = sum(y)/P.Vres;
    y  = y/Sy;
end


% plot(y)% plot([y' y2'])
% plot([F' G'])
end