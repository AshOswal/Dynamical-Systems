function simulate_nonlinear (input_E, input_I)

% function simulates a circuit composed from excitatory and inhibitory
% populations without delay but with self-excitatory feedback
% potentialy receiving slow excitatory inputs
% Optional inputs:
%  input_E - flag indicating if excitatory population gets slow input
%            (default 1)
%  input_I - flag indicating if inhibitory population gets slow input
%            (default 0)

if nargin < 1
    input_E = 0.05;
end
if nargin < 2
    input_I = 0;
end
%% constants
DT = 0.01;      %integration step
MAXT = 0.2;     %duration of simulation
wI2E = 2;       %weight of connections from inhibitory to excitatory
wE2I = 2;       %weight of connections from excitatory to inhibitory
wE2E = 2;%3;       %weight of connection from excitatory to excitatory; 2 is the value used by Angela, but 3 necessary for PAC
InputE = 0.5;    %constant input to excitatory population; value of around 0.5 required for bringing into linear range
InputI = 0.3;
%freq = 0.020;      %frequency of slow oscillations
freq = 25;
%% simulation
% initial condition
E(1) = InputE;
I(1) = InputI;

%{
for index = 1:MAXT/DT-1
    %integration of excitatory population
    Total_input_E(index) = InputE;
    if input_E
        Total_input_E(index) = Total_input_E(index) ;%+ (sin(index*DT*2*pi*freq)+1)/4;
    end    
    
    % Euler's method for integration
    E(index+1) = E(index) + DT * (sigmoid(Total_input_E(index) + wE2E*E(index) - wI2E*I(index)) - E(index));
    if E(index+1) < 0
        E(index+1) = 0;
    end
    
    %integration of inhibitory population
    Total_input_I(index) = InputI;
    if input_I
        Total_input_I(index) = Total_input_I(index) + (sin(index*DT*2*pi*freq)+1)/4;
    end    
    
    % Euler's method for integration
    I(index+1) = I(index) + DT * (sigmoid(Total_input_I(index) + wE2I*E(index)) - I(index));
    if I(index+1) < 0
        I(index+1) = 0;
    end
end
%}
%% phase portrait

% y(1) = E;
% y(2) = I;
c_e = 0.5;
c_i = 0;
i_e = @(t) InputE + c_e*(sin(t*2*pi*freq)+1);
i_i = @(t) InputI + c_i*(sin(t*2*pi*freq)+1)/4;

step = 0.1;
f = @(t,y)[-y(1)+sigmoid(i_e(t) + (wE2E*(y(1))) - (wI2E*(y(2))));-y(2)+sigmoid(i_i(t) + wE2I*y(1))]./(0.0032/4);

[xx yy] = meshgrid(0:step:1,0:step:1);

u = zeros(size(xx));
v = zeros(size(xx));

% we can use a single loop over each element to compute the derivatives at
% each point (y1, y2)
t=0; % we want the derivatives at each point at t=0, i.e. the starting time
for i = 1:numel(xx)
    Yprime = f(t,[xx(i); yy(i)]);
    Yprime = Yprime./norm(Yprime);
    u(i) = Yprime(1); % E
    v(i) = Yprime(2); % I
end

% nullcline functions - determine I as a function of E
Ev = 0:step/100:1;
Iv = 0:step/100:1;
Enull =  (E(1) + wE2E.*Ev - 1 + (0.25*log((-Ev+1)./Ev)))/wI2E;
Inull = 1./(1 + exp(-4*I(1) - 4*wE2I.*Ev +4));

% integrate using ode45 for a fixed time window
[time, ei]  = ode45(f,[0 MAXT],[0 0]);
figure;quiver(xx,yy,u,v,'k');hold on;
plot(Ev,Inull,'r-');hold on
plot(Ev,Enull,'b-');hold on
plot(ei(:,1),ei(:,2),'y');
xlabel('Excitatory')
ylabel('Inhibitory')
axis tight equal;


figure;
subplot(2,1,1)
plot(time,[i_e(time),i_i(time)]);
title ('Input to the populations');
xlabel ('Time');
legend ('Excitatory', 'Inhibitory');

subplot(2,1,2);
plot(time,ei);
title ('Activity of the populations');
xlabel ('Time');
legend ('Excitatory', 'Inhibitory');


%% sigmoid function and its derivative

y = sigmoid(0:0.1:3);
ydiff = 4.*y.*(1-y);
%figure;plot(0:0.1:3,y);hold on;plot(0:0.1:3,ydiff,'r');

%{
%% visualization
figure;
subplot(2,1,1)
plot (DT:DT:MAXT-DT, Total_input_E);
hold on
plot (DT:DT:MAXT-DT, Total_input_I, 'r');
title ('Inputs to the populations');
xlabel ('Time');
legend ('Excitatory', 'Inhibitory');

subplot(2,1,2);
plot (DT:DT:MAXT,E);
hold on
plot (DT:DT:MAXT,I, 'r');
title ('Activity of the populations');
xlabel ('Time');
legend ('Excitatory', 'Inhibitory');
%}
end

%% sigmoid function
function y = sigmoid (x)
y = 1./(1+exp(-4*(x-1)));
end

