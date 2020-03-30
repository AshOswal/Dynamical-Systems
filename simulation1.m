function simulation1



% simulation 1 discussed with Rafal

p = [];
p.dt          =  0.01;  % time step [ms] 
p.t_Start     =  0;     % simulation start time [ms] 
p.t_End       =  1000;   % simulation end time [ms]
p.E_L         = -70;    % resting membrane potential [mV]
p.V_th        = -55;    % spike threshold [mV]
p.V_reset     = -75;    % value to reset voltage to after a spike [mV] 
p.V_spike     =  20;    % the potential of the spike
p.R_m         =  10;    % the membrane resistance %
p.tau         =  10;    % membrane time constant [ms]
p.I_threshold = (p.V_th - p.E_L)/p.R_m; % current below which cell does not fire
p.Gc          = p.I_threshold/100;     
p.popN        = 100;    % Number of IF neurones in population model 
p.fGPall      = 0:0.025:0.1; % Firing rate [Hz]*10-3 
p.fCtx        = 2;      % Firing rate [Hz]*10-3
p.dGP         = 0;      % delay term [ms]
p.dCtx        = 0;      % delay term [ms]



V             = zeros(p.popN,numel(p.t_Start:p.dt:p.t_End),numel(p.fGPall));
V(:,1,:)      = p.E_L;
Vv            = V;
spike         = zeros(p.popN,numel(p.fGPall));
    
for f = 1:numel(p.fGPall)
   
    p.fGP         = p.fGPall(f);
    GP            = sum(simulate_poisson_train(p.popN,p.fGP,p.t_Start,p.t_End,p.dt),1);     % Inhibitory
    Ctx           = sum(simulate_poisson_train(p.popN,p.fCtx,p.t_Start,p.t_End,p.dt),1);    % Excitatory
    i_GP          = -GP*p.Gc/p.dt;  
    i_Ctx         = Ctx*p.Gc/p.dt; 

    for N = 1:p.popN

        for t = 2:size(V,2)

            % integral of IF neurone
            Ie       = i_GP(t-1) + i_Ctx(t-1);

            V(N,t,f) = p.E_L + p.R_m*Ie + (V(N,t-1,f) - p.E_L - p.R_m*Ie)*exp(-p.dt/p.tau);

            % reset and make spike

            if  V(N,t,f)   > p.V_th
                V(N,t,f)   = p.V_reset;
                Vv(N,t,f)  = p.V_spike;
                spike(N,f) = spike(N,f) + 1;         
            else
                Vv(N,t,f) = V(N,t,f);
            end

        end

        fHz(N,f) = spike(N,f)/(p.t_End - p.t_Start);

    end

end

figure;plot(p.fGPall,mean(fHz,1),'ro');
% predicted and actual slopes and intercepts

i = find(mean(fHz,1) > 0); 
obs_b  = regress(mean(fHz(:,i),1)',p.fGPall(i)');
pred_b = (p.R_m*p.Gc*p.popN)/(p.tau*(p.V_th-p.V_reset));

obs_i  = mean(fHz(:,1),1);
pred_i =  (p.R_m*p.Gc*p.popN*p.fCtx)/(p.tau*(p.V_th-p.V_reset)) + ((p.E_L - p.V_th)/(p.tau*(p.V_th-p.V_reset))); %Ignore this term;

title(sprintf('predicted slope %0.3f \n actual slope %0.3f \n predicted intercept %0.3f \n actual intercept %0.3f',pred_b,obs_b,pred_i,obs_i));




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
