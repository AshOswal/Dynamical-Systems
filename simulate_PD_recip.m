function simulate_PD_recip
clc;
close all;
clear all;

% the parameter descriptions of the original model as per Nevado Holgado 2010
D = 1;
p.dTsg = D*6e-3;   % ms
p.dTgs = D*6e-3;   % ms
p.dTgg = D*4e-3;   % ms
p.dTcs = D*5.5e-3;

p.Ts   = 6e-3;   % ms
p.Tg   = 14e-3;  % ms

p.Ctx  = 27;     % spm/s
p.Str  = 2;      % sp/s
p.Ms   = 300;    % sp/s
p.Bs   = 17;     % sp/s
p.Mg   = 400;    % sp/s
p.Bg   = 75;     % sp/s

% for the cortical populations
p.K = 1;
% disease
p.Wsg  = 19 + p.K*(20-19);
p.Wgs  = 1.12 + p.K*(10.7-1.12);
p.Wgg  = 6.6 + p.K*(12.3 - 6.6);
p.Wcs  = 2.42 + p.K*(9.2 - 2.42);
p.Wxg  = 15.1 + p.K*(139.4 - 15.1);

% plot sigmoid relationships
Fs = sigmoid(-200:700,p.Ms,p.Bs);
Fg = sigmoid(-200:700,p.Mg,p.Bg);
figure;subplot(2,1,1);plot(Fs);subplot(2,1,2);plot(Fg);

% integrate model
IC =    zeros(2,1);
ops_DDE         = ddeset('RelTol',1e-6,'AbsTol',1e-6);
sol             = dde23(@ddefun, [p.dTsg,p.dTgs,p.dTgg]', IC, [0 1], ops_DDE,p);
sol.y           = sol.y; %+ randn(size(sol.y,1),size(sol.y,2));

figure;%subplot(2,1,1);
plot(sol.x,sol.y(1,:),'b');hold on;
plot(sol.x,sol.y(2,:),'g');hold on;
legend({'STN','GP','Cortex E','Cortex I'});
xlabel('time');ylabel('rate');
title('DDE 23');
end


function y = sigmoid (x,max,min)
y = max./(1+(exp(-4*x*inv(max))*(max-min)/min));
end

function dydt = ddefun(t,Y,Z,p)
% the integral function
% Y(1) = STN
% Y(2) = GP

% Z(:,1) = dTsg
% Z(:,2) = dTgs
% Z(:,3) = dTgg
% Z(:,4) = dTcs

% remember that Z is coded Z(i,j) such that i is the giving population and
% j is the delay

dydt    = zeros(2,1);
dydt(1) =  (sigmoid((-p.Wgs*(Z(2,2)) + p.Wcs*p.Ctx),p.Ms,p.Bs) - Y(1))/p.Ts; % STN
dydt(2) =  (sigmoid(( p.Wsg*(Z(1,1)) - p.Wgg*(Z(2,1)) - p.Wxg*p.Str),p.Mg,p.Bg) - Y(2))./p.Tg; % GP

end

