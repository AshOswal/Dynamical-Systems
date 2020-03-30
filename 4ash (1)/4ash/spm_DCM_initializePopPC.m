function [Ep] = spm_DCM_initializePopPC(M0,pop,Y)
%% Test Parameters - From BeuroElectro.com
% TestParam7
% clear all
%% Priors
% --------------------
% --- Expectations ---
% --------------------

P.ang    = M0.P.ang;

% --- Conectivity matrices ---

P.A(1).M = [ -32;
             -32];
    
P.A(2).M = [-32;
            -32];
        
P.A(3).M = [ -32;
             -32];
    
P.A(4).M = [-32;
            -32];
        
P.C(1).M = [log(.1)]';

P.C(2).M = [log(.1)]';

P.C(3).M = [-32]';

P.C(4).M = [-32]';

% --- inhibitory interneurons ---
P.P(1) = M0.P.P(pop);

CP.ang    = 10*pi;

CP.A(1).M = exp([ -32;
                  -32]);
    
CP.A(2).M = exp([-32;
                 -32]);

CP.A(3).M = exp([ -32;
                  -32]);
    
CP.A(4).M = exp([ -32;
             	  -32]);

CP.C(1).M = exp([4]');
CP.C(2).M = exp([4]');

CP.C(3).M = exp([-32]');

CP.C(4).M = exp([-32]');

% --- inhibitory interneurons ---

CP.P(1) = M0.pC.P(pop);


%%
M.opt.N              = 4000;
M.opt.r              = 1;
M.opt.popCurrents    = [2 10; 2 -70; 1 10; 1 -70];
M.opt.popFiringRates = [2; 1];
M.opt.dpsize         = 10^-6;
M.opt.svdtol         = 10^-6;
M.opt.Drawplots      = 0;
M.opt.vmin           = -16;
M.opt.vmax           = 8;
M.opt.Nmax           = 32;
M.opt.sigma          = 0; 
M.opt.Ncoefs = length(Y.y);
Vp = 0*spm_vec(P); Vp = [Vp Vp]; Vp(4,1) = 1; Vp(4,2) = -1;
% Vp     = spm_svd(diag(spm_vec(CP)),exp(-32));
%% Data 

M.pE = P;
M.pC = CP;
M.hE = 0*[3 9 -3 -3 3 3 3 9]';
M.hC = eye(length(M.hE))*.01;



YGI =  Y.y(:,1)/(M.opt.popCurrents(1,2) - P.P(1).Vg(2));
YGE =  Y.y(:,2)/(M.opt.popCurrents(2,2) - P.P(1).Vg(1));
YGI = YGI - min(YGI);
YGE = YGE - min(YGE);

% figure
% plot([YGE YGI])

YGIf = fft(YGI);
YGEf = fft(YGE);
Coefs = [YGEf'; YGIf';zeros(2, length(YGI))];

f0  = Y.y(1,5);
N = length(Coefs);
Vfreq = mod((0:(N-1))+N/2,N)-N/2;

Uf = @(t) 1/N*real(Coefs*exp(Vfreq'*(-1i*2*pi*f0*t)));
figure
vt = 1/N*1/f0:1/N*1/f0/4:1/f0;
vt0  = 1/N*1/f0:1/N*1/f0/1:1/f0;
plot(vt,Uf(vt)); hold on
plot(vt0, [YGE YGI],'-')

M.opt.U.f0 = Y.y(1,5);
M.opt.U.u = Uf;

M.opt.popCurrents = [];
M.opt.popFiringRates = 1;
% -------------------------------------------
popfiringrate = 6;
M.hE = M.hE(popfiringrate)+0;
M.hC = M.hC(popfiringrate,popfiringrate);
Y.y  = Y.y(:,popfiringrate); % First population firing rate
% -------------------------------------------
M.opt.Drawplots=0;
M.converged = 0;
M.F = -Inf;
M.opt.LCtol = .1;
M.opt.Drawprogress = 1;
% M.opt.svdtol = 0;
% M.opt.dpsize       = 10^-6;
M.opt.vmin         = -16;
M.opt.vmax         = 10;
M.opt.vini         = -4;
M.opt.T0max        = 1.5/f0;


[Ep,Cp,Eh,F,dFdp,dFdpp,f] = spm_nlsi_GN_LC_fr_spinit(M,0,Y);

end





