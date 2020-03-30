clear all;

% number of sources and LFP channels (usually the same)
%--------------------------------------------------------------------------
n     = 1; % number of sources
nc    = 1; % number of channels

% specify model
%--------------------------------------------------------------------------
M.f   = 'spm_fx_bg_s';
M.x   = sparse(n,2);
M.Nit = 2^8;
%M.D = 1;
%M.n  = n*13;
%M.pE = pE;
%M.pC = pC;
%M.m  = n;
%M.l  = nc;
%M.Hz = (1:64)';
 
P = [];

% integrate with no inputs 
%--------------------------------------------------------------------------
N    = 100;
U.dt = 1/100;
U.u  = zeros(1,N)';
%U.u  = randn(N,M.m)/16;
%U.u  = sqrt(spm_Q(1/16,N))*U.u;
LFP  = spm_int_L(P,M,U);
figure;
plot((1:N).*U.dt,LFP(:,1),'b');hold on;
plot((1:N).*U.dt,LFP(:,2),'g');hold on;
%plot((1:N).*U.dt,LFP(:,3),'r');hold on;
%plot((1:N).*U.dt,LFP(:,4),'k');hold on;
xlabel('time');
ylabel('firing rate');
legend({'STN','GP','E','I'});
title('SPM INT L');