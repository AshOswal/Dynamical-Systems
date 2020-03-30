function simulateICP
P = [];
clear all;
close all;

% state variables which will change over time dfdt
P.Pic   = 9.5;        % mmHg ICP
P.Ca    = 0.15;       % ml/mmHg Compliance

% input variables
a              = repmat(100,1,2000);
b              = repmat(100,1,2500);
c              = repmat(100,1,10*1e3);
P.Pa           = [a c b];                                                                                   % mmHg Arterial Pressure
P.T            = [1:length(a), length(a)+1e-3:1e-3:length(a)+10, length(a)+11:length(a)+10+length(b)];      % small temporal step size for integration

P.dPadt        = [0 diff(P.Pa)];

    

P.vs    = 6.0;        % mmHg dural venous pressure
P.I     = 0;          % ml/s CSF infusion rate

% constant parameters
P.Ro    = 6*526.3;      % mmHg x s/ml CSF outflow resistance
P.Rpv   = 1.24;       % mmHg x s/ml Proximal venous resistance
P.Rf    = 2.38*1e3;   % mmHg s s/ml CSF fluid formation resistance
P.DCa1  = 0.75;       % ml/mmHg Sigmoid curve amplitude
P.DCa2  = 0.075;      % ml/mmHg Sigmoid curve amplitude 
P.Can   = 0.15;       % ADD HERE
P.Kr    = 4.91*1e4;   % mmHg3 x s /ml;
P.KE    = [ 2.1 0.44 0.27];       % 
P.t     = 20;         % seconds
P.Qn    = 12.5;       % ml/s
P.G     = 1.5;        % ml x /mmHg x 100% /CBF change
P.Pc    = 25;
Z       = {P.Pa P.dPadt P.vs P.I};

figure;
col = {'k','b','r','g','m','c'};
for j = 1:length(P.KE)
    P.Ke = P.KE(j);
    [T, f]  = ode15s(@(t,F) odefx(t,F,Z,P),P.T,[P.Pic P.Ca]);
    VolArt  = f(:,2)'.*(P.Pa - f(:,1)');
    Ra      = [];
    
    for t = 1:length(T)
        [~,r] = odefx(T(t),f(t,:),Z,P);
        Ra    = [Ra r];
    end
    CBF     = (P.Pa - P.Pc)./Ra;
    
    subplot(3,1,1)
    plot(T,f(:,1),col{j});xlabel('time (s)');ylabel('ICP mmHg');hold on;
    legend({'ICP'});
    %axis([0 300 0 50]);
    
    subplot(3,1,2)
    plot(T,CBF,col{j});xlabel('time (s)');ylabel('CBF ml/s');hold on;
    legend('CBF');
    %axis([0 300 0 30]);
    
    subplot(3,1,3)
    plot(VolArt,f(:,1),col{j});xlabel('Art Vol');ylabel('ICP');hold on;
    legend('CBF');
    
end

function [dfdt,Ra] = odefx(t,F,Z,P)
% initiate equations of motion

% Z(1)    = P.Pa;
% Z(2)    = P.dPadt;
% Z(3)    = P.vs;
% Z(4)    = P.I;

% F(1)    = P.Pic;
% F(2)    = P.Ca;

% dfdt(1) = dCadt;
% dfdt(2) = dPicdt;

z1      = interp1(P.T,Z{1},t);
z2      = interp1(P.T,Z{2},t);


Va      = F(2)*(z1 - F(1));
Ra      = P.Kr * (P.Can^2)/(Va^2);
Pc      = ((z1*P.Rpv) + (F(1)*Ra))/(P.Rpv +Ra);
q       = (z1 - Pc)/Ra;
x       = (q - P.Qn)/P.Qn;
if x<0
   DCa  = P.DCa1;   
else
   DCa  = P.DCa2;
end
Ksig    = DCa/4;

SigGx   = ((P.Can + DCa/2) + (P.Can - DCa/2)*exp(P.G*x/Ksig)) / (1+ exp(P.G*x/Ksig));
%SigGx   = Sig(P,x,Ksig,DCa);


dCadt   = (1/P.t)*(-F(2) + SigGx);

t1 = (P.Ke*F(1))/(1 + (F(2)*P.Ke*F(1)));
t2 = F(2)*z2;
t3 = dCadt*(z1 - F(1));
t4 = (Pc - F(1))/P.Rf;
t5 = (F(1)- Z{3})/P.Ro;

dPicdt  = t1*(t2+t3+t4-t5+Z{4});
dfdt    = [dPicdt;dCadt];



function s = Sig(P,x,Ksig,DCa)

s = ((P.Can + DCa/2) + (P.Can - DCa/2)*exp(P.G*x/Ksig)) / (1+ exp(P.G*x/Ksig));
