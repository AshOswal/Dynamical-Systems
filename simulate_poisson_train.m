function spike = simulate_poisson_train(N,r,ts,tf,dt)
% use simple approach suggested by Dayan and Abbott

% N - number of neurones
% r - rate of poisson firing 1/isi [Hz - 1/s] 
% ts - start time                  [s]
% tf - finish time                 [s]
% dt - time step

out = rand(N,((tf-ts)/dt)+1);
spike = (r*dt)> out;
