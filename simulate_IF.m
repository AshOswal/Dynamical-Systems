function simulate_IF

% simulate IF neuron receiving spike inputs from 100 neurones. For
% simplicity in this case we assume that each input multiplied by the
% corresponding synaptic conductance is suprathreshold. This avoids having
% to consider what happens to multiple subthreshold inputs occuring
% simultaneously - i.e maybe they are additive - we assume that they are
% for suprathreshold inputs


%% DEFINE NEURONAL PARAMETERS
clear all;
dt          =  0.01;  % time step [ms]
t_Start     =  0;     % simulation start time [ms] 
t_End       =  100;   % simulation end time [ms]
E_L         = -70;    % resting membrane potential [mV]
V_th        = -55;    % spike threshold [mV]
V_reset     = -75;    % value to reset voltage to after a spike [mV] 
V_spike     =  20;    % the potential of the spike
R_m         =  10;    % the membrane resistance
tau         =  10;    % membrane time constant [ms]
I_threshold = (V_th - E_L)/R_m; % current below which cell does not fire
Gc          = I_threshold;      
popN        = 100;    % Number of IF neurones in population model 

frate = 1:10; % firing rates to simulate ms-1
truerate = (1000/(1/dt)).*frate;
AveRate = nan(popN,numel(frate));
%%

t_vect = t_Start:dt:t_End;                                    % will hold vector of times
V_vect = repmat([E_L zeros(1,length(t_vect)-1)],popN,1);
V_plot_vect = repmat([E_L zeros(1,length(t_vect)-1)],popN,1); % plotting that displays a spike

for Fr  = frate; % loop over different input firing rates
    
    input_spike_train = simulate_poisson_train(1000,Fr,t_Start,t_End,dt);

    for n = 1:popN
        
    jj  =     input_spike_train(randperm(100),:);
    I_e_vect = sum(jj,1).*Gc;     % sum sparse connectivity and multiply by conductance    
    
    i = 1;
    Nspike = 0;
    for t = dt:dt:t_End
        V_vect(n,i+1) = E_L + R_m*I_e_vect(i) + (V_vect(n,i) - E_L - R_m*I_e_vect(i))*exp(-dt/tau);   % analytical integral of IF      
        % detect and reset spike
        if  V_vect(n,i+1) > V_th
            V_vect(n,i+1) = V_reset;
            V_plot_vect(n,i+1) = V_spike;
            Nspike = Nspike + 1;
        else
            V_plot_vect(n,i+1) = V_vect(n,i+1);
        end       
        i = i+1;
    end
    
    AveRate(n,find(Fr==frate)) = Nspike/(t_End - t_Start);
    %{
    figure(1);
    subplot(length(frate),1,Nplot);
    plot(t_vect, V_plot_vect);
    title('Voltage vs. time');
    xlabel('Time in ms');
    ylabel('Voltage in mV');
    %}
    end
end

% we need to compare actual and predicted firing rate
figure;
plot(truerate,mean(AveRate,1),'b--')
xlabel('Firing rate of input neurones (Hz)')
ylabel('Firing rate of output neurone (Hz)')

% compare actual and predicted gradients since these corespond to the
% weights
wpred = (R_m*Gc*popN)/(tau*(V_th-V_reset));
wactual = (max(mean(AveRate,1))-min(mean(AveRate,1)))/(max(truerate)-min(truerate));
title(sprintf('predicted firing rate %0.2f \n actual firing rate %0.2f',wpred,wactual));


function spike = simulate_poisson_train(N,r,ts,tf,dt)
% use simple approach suggested by Dayan and Abbott

% N - number of neurones
% r - rate of poisson firing 1/isi [Hz]
% ts - start time
% tf - finish time
% dt - time step

out = rand(N,((tf-ts)/dt)+1);
spike = (r*dt)> out;
% the homogeneous poisson spike generator

