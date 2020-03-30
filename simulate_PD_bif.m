function [Y,F,S] = simulate_PD_bif(varargin)
clc;
close all;
%clear all;

% the parameter descriptions of the original model as per Nevado Holgado 2010

%%

% note we add autoinhibition in cortical layers in accordance with the CMC
% paper of Rosalyn Moran. Additional additions are a model of the GPi. The
% time constants of the GPi are taken from the literature. 


if isempty(varargin)
    p.D        = 1;
    D = p.D;
    p.dTsg     = 6e-3;%D*6e-3;   % ms
    p.dTgs     = 6e-3;%D*6e-3;   % ms
    p.dTgg     = D*4e-3;   % ms
    p.dTcs     = D*5.5e-3;
    p.dTsc     = D*20e-3;  %21.5e-3;
    p.dTsgpi   = 6e-3;
    p.dTgpegpi = 6e-3;
    p.dTgpigpi = 4*1e-3;
    
    p.Ts   = 12.8*1e-3;   % ms %12.8
    p.Tg   = 20*1e-3;     % ms %20
    p.Te   = 11e-3;%10e-3;  % ms
    p.Ti   = 11e-3;%10e-3;  % ms
    p.Tgpi = 14e-3;       % see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2585404/#!po=67.1053 for reference 
    
    p.Ctx  = 172.18; %27;       % spm/s
    p.Str  = 8.46;   %  2;      % sp/s
    p.Ms   = 300;               % sp/s
    p.Bs   = 10;%17;     % sp/s
    p.Mg   = 400;    % sp/s
    p.Bg   = 20;%75;     % sp/s
    p.Mgi  = 300;
    p.Bgi  = 18;
    p.Mi   = 205;    % normal range 200-330 sp/s
    p.Bi   = 9.87;     % normal range 0-20 sp/s
    p.Me   = 75;     % normal range 50-80 sp/s
    p.Be   = 17.85;     % normal range 0-20 sp/s
    
    % for the cortical populations
    p.dTcs  = D*5.5e-3;   % ms - allowed range 1-10 ms
    p.dTcc  = D*6e-3; %5e-3%1e-3;%4.65e-3;  % 2;     % ms - allowed range 1-10 ms
    
    p.K = 1;
    % disease
    p.Wsg  = 4.87;%19 + p.K*(20-19);
    p.Wgs  = 1.33;%1.12 + p.K*(10.7-1.12);
    p.Wgg  = 0.53;%0.53;%6.6 + p.K*(12.3 - 6.6);
    p.Wcs  = 10;%2.42 + p.K*(9.2 - 2.42);
    
    % cortical parameters
    p.Wsc     = 0;%0.4;%0.4;%0.4;%0.4%8.93;
    p.Wcc     = 4;
    p.Wsgpi   = 3.7;
    p.Wgpegpi = 1;
    
    % stochastic component - GWN
    p.S       = 0;
    
else
    p = varargin{1};
end

% integrate model
IC              =    zeros(2,1);
ops_DDE         = ddeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',5*1e-3);
%ops_DDE         = ddeset('MaxStep',1e-3);
sol             = dde23(@ddefun, [p.dTsg,p.dTgs,p.dTgg]', IC, [0 3], ops_DDE,p);
%sol             = dde23(@ddefun, [p.dTsg,p.dTgs,p.dTgg, p.dTcs, p.dTsc, p.dTcc]', IC, [0 2], ops_DDE,p);

sol.y           = sol.y + 0*randn(size(sol.y,1),size(sol.y,2));

figure;subplot(3,1,1);
plot(sol.x,sol.y(1,:),'b','LineWidth',1);hold on;
plot(sol.x,sol.y(2,:),'g','LineWidth',1);hold on;
%plot(sol.x,sol.y(3,:),'r','LineWidth',1);hold on;
%plot(sol.x,sol.y(4,:),'k','LineWidth',1);hold on;
%plot(sol.x,sol.y(5,:),'m','LineWidth',1);hold on;
axis([0 3 0 max(sol.y(:))+20]);
%h = legend({'STN','GPe','GPi','Cortex E','Cortex I'});
%set(h,'box','off','LineWidth',2);
xlabel('Time (s)');ylabel('Firing Rate spk/s');
set(gca,'FontSize',12);
%title('DDE 23');


[S,n,f] = ft_specest_mtmfft(sol.y(:,sol.x>1 & sol.x<4),sol.x(sol.x>1 & sol.x <4),'taper','dpss','tapsmofrq',1,'freqoi',[0:1: 50]);
s = squeeze(mean(abs(S),1));
subplot(3,1,2);
plot(f,s(1,:),'b','LineWidth',2);hold on;
plot(f,s(2,:),'g','LineWidth',2);hold on;
%plot(f,s(3,:),'r','LineWidth',2);hold on;
%plot(f,s(4,:),'k','LineWidth',2);hold on;
%plot(f,s(5,:),'m','LineWidth',2);hold on,
h = legend({'STN','GPe'});
set(h,'box','off','LineWidth',2);
xlabel('Frequency (Hz)');ylabel('Power (AU)');
if isfield(p,'axis')
axis([0 p.axis(1) 0 p.axis(2)]);
end
set(gca,'FontSize',12);

subplot(3,1,3)
plot(sol.y(1,:),sol.y(2,:),'k')
xlabel('STN rate (Hz)');
ylabel('GPe rate (Hz)');
box on

Y.y = sol.y;
Y.x = sol.x;
F.f = f;
F.s = s;


end


function y = sigmoid (x,max,min)
y = max./(1+(exp(-4*x*inv(max))*(max-min)/min));
end

function dydt = ddefun(t,Y,Z,p)
% the integral function
% Y(1) = STN
% Y(2) = GP
% Y(3) = GPi
% Y(4) = Cortex - E
% Y(5) = Cortex - I

% Z(:,1) = dTsg
% Z(:,2) = dTgs
% Z(:,3) = dTgg
% Z(:,4) = dTcs
% Z(:,5) = dTsc  % Note that this is now an inhbitory effect of the GPi
% Z(:,6) = dTcc
% Z(:,7) = dTsGPi
% Z(:,8) = dTGPeGPi
%

% remember that Z is coded Z(i,j) such that i is the giving population and
% j is the delay

dydt = zeros(2,1);


% dydt(1) =   (sigmoid(-p.Wgs*Z(2,2)   + p.Wcs*Z(3,4),p.Ms,p.Bs) - Y(1)) /p.Ts  ;      % STN
% dydt(2) =   (sigmoid( p.Wsg*Z(1,1)   - p.Wgg*Z(2,3) - p.Str,p.Mg,p.Bg) - Y(2))/p.Tg; % GP
% dydt(3) =   (sigmoid(-p.Wsc*Z(1,5)   - p.Wcc*Z(4,6) + p.Ctx,p.Me,p.Be) - Y(3))/p.Te; % Cortex E
% dydt(4) =   (sigmoid(p.Wcc*Z(3,6)    - p.Wcc*Z(4,6),p.Mi,p.Bi) - Y(4))/p.Ti;         % Cortex I
S       =   p.S;
%dydt(1) =   (sigmoid(-p.Wgs*Z(2,2)   + p.Wcs*Z(4,4),p.Ms,p.Bs)                                 - Y(1) + S*randn)/p.Ts  ; % STN
dydt(1) =   (sigmoid(-p.Wgs*Z(2,2)   + p.sinput,p.Ms,p.Bs)                                 - Y(1) + S*randn)/p.Ts  ; % STN

dydt(2) =   (sigmoid( p.Wsg*Z(1,1)   - p.Wgg*Z(2,3) - p.Str,p.Mg,p.Bg)                         - Y(2) + S*randn)/p.Tg;   % GPe
%dydt(3) =   (sigmoid( p.Wsgpi*Z(1,6) - p.Wgpegpi*Z(2,7) - p.Wgg*Z(3,8) -p.Str,p.Mgi,p.Bgi)       - Y(3) + S*randn)/p.Tgpi; % GPi          
%dydt(4) =   (sigmoid(-p.Wsc*Z(3,4)   - p.Wcc*Z(5,5) + p.Ctx,p.Me,p.Be)                         - Y(4) + S*randn)/p.Te; % Cortex E
%dydt(5) =   (sigmoid(p.Wcc*Z(4,5)    - p.Wcc/1*Z(5,5),p.Mi,p.Bi)                                 - Y(5) + S*randn)/p.Ti;         % Cortex I

end

