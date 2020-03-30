function simulate_PD_betaPAC_IF
clear all;
close all;
p = [];
p.dt          =  0.01;  % time step [ms]
p.t_Start     =  0;     % simulation start time [ms]
p.t_End       =  1000;  % simulation end time [ms]
p.E_L         = -70;    % resting membrane potential [mV]
p.V_th        = -55;    % spike threshold [mV]
p.V_reset     = -75;    % value to reset voltage to after a spike [mV]
p.V_spike     =  20;    % the potential of the spike
p.R_m         =  10;    % the membrane resistance %
p.tau         =  10;    % membrane time constant [ms]
p.I_threshold = (p.V_th - p.E_L)/p.R_m; % current below which cell does not fire
p.Gc          = p.I_threshold; % may need diferent terns for differnt synapses
p.popN        = 500;    % Number of IF neurones in population model
p.popC        = 100;    % Number of inputs received by each neurone - sparse connectivity
p.fStr        = 2*1e-3; % Firing rate [Hz]*10-3  -  see Nevado Holgado et al. 2010
p.fCtx        = 27*1e-3;% Firing rate [Hz]*10-3  -  see Nevado Holgado et al. 2010
p.tGP         = 8;      % delay term [ms]        -  see Nevado Holgado et al. 2010
p.tSTN        = 8;      % delay term [ms]        -  see Nevado Holgado et al. 2010

V_STN         = zeros(p.popN,numel(p.t_Start:p.dt:p.t_End));
V_GP          = zeros(p.popN,numel(p.t_Start:p.dt:p.t_End));

V_STN(:,1)    = p.E_L;
V_GP(:,1)     = p.E_L;
Vv_STN        = V_STN;
Vv_GP         = V_GP;

Str_all           = (simulate_poisson_train(p.popN,p.fStr,p.t_Start,p.t_End,p.dt));     % Inhibitory
Ctx_all           = (simulate_poisson_train(p.popN,p.fCtx,p.t_Start,p.t_End,p.dt));     % Excitatory
i_Str_all         = -Str_all*p.Gc/p.dt;
i_Ctx_all         = Ctx_all*3*p.Gc/p.dt;

STN_spike         = zeros(p.popN,numel(p.t_Start:p.dt:p.t_End));
GP_spike          = zeros(p.popN,numel(p.t_Start:p.dt:p.t_End));

rnd_mat           = [];
for n = 1:p.popN
    rnd_mat = [rnd_mat,randperm(p.popN,p.popC)']; % code for generating fixed pattern of sparse connectivity
end

for t = 2:size(V_STN,2)
    
    for N = 1:p.popN % must integrate populations before time steps
        
        % integrate STN inputs assuming *sparse* connectivity from inputs and EI neurones as per Brunel 1999     
        delay = p.tGP*1/p.dt;
        if t > delay
            I_STN      = -sum(GP_spike(rnd_mat(:,N),t-delay)*p.Gc/p.dt,1) + sum(i_Ctx_all(rnd_mat(:,N),t-1),1);
        else
            I_STN      = sum(i_Ctx_all(rnd_mat(:,N),t-1),1);
        end
        
        V_STN(N,t) = p.E_L + p.R_m*I_STN + (V_STN(N,t-1) - p.E_L - p.R_m*I_STN)*exp(-p.dt/p.tau);
        
        % reset and make spike
        if  V_STN(N,t)     > p.V_th
            V_STN(N,t)     = p.V_reset;
            Vv_STN(N,t)    = p.V_spike;
            STN_spike(N,t) = 1;
        else
            Vv_STN(N,t) = V_STN(N,t);
        end
        
        % integrate GP inputs
        delay = p.tSTN*1/p.dt;
        if t > delay
            I_GP      = sum(STN_spike(rnd_mat(:,N),t-delay)*0.5*p.Gc/p.dt,1) - sum(i_Str_all(rnd_mat(:,N),t-1),1);
        else
            I_GP      = -sum(i_Str_all(rnd_mat(:,N),t-1),1);
        end
        
        V_GP(N,t) = p.E_L + p.R_m*I_GP +  (V_GP(N,t-1) - p.E_L - p.R_m*I_GP)*exp(-p.dt/p.tau);
        
        % reset and make spike
        if  V_GP(N,t)     > p.V_th
            V_GP(N,t)     = p.V_reset;
            Vv_GP(N,t)    = p.V_spike;
            GP_spike(N,t) = 1;
        else
            Vv_GP(N,t) = V_GP(N,t);
        end

    end
   
    f_STN(N) = sum(STN_spike(N,:))/(p.t_End - p.t_Start);
    f_GP(N)  = sum(GP_spike(N,:))/(p.t_End - p.t_Start);
end

% compute firing rates over a gussian window for smoothness as per Dayan
% and Abbott
sigma = 200;
width = round((6*sigma - 1)/2);
gauss = normpdf(-width:width,0,sigma);
gauss = gauss./sum(gauss);

for n = 1:p.popN
    STN_spike_smooth(n,:) = conv(STN_spike(n,:),gauss,'same');
    GP_spike_smooth(n,:)  = conv(GP_spike(n,:),gauss,'same');
end

%% response of single neurones in STN and GP
figure;
subplot(2,2,1);
plot(p.t_Start:p.dt:p.t_End,STN_spike');
title('STN spiking');

subplot(2,2,2);
plot(p.t_Start:p.dt:p.t_End,GP_spike');
title('GP spiking');

subplot(2,2,3);
plot(p.t_Start:p.dt:p.t_End,STN_spike_smooth');
title('STN firing rate');

subplot(2,2,4);
plot(p.t_Start:p.dt:p.t_End,GP_spike_smooth');
title('GP firing rate');



function spike = simulate_poisson_train(N,r,ts,tf,dt)

% use simple approach suggested by Dayan and Abbott for a homogeneous
% poisson process

% N - number of neurones
% r - rate of poisson firing 1/isi
% ts - start time
% tf - finish time
% dt - time step
out = rand(N,((tf-ts)/dt)+1);
spike = (r*dt)> out;

