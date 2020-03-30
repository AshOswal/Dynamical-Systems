function test

a =0.25;
b= -0.01;
c= -1.00;
d=0.01;

y10 = 80; % initial population of prey
y20 = 30; % initial population of predator
tspan = [0 200]; % time span to integrate over
[t,Y] = ode45(@myode,tspan,[y10 y20]);

y1 = Y(:,1);
y2 = Y(:,2);

figure;
plot(y1,y2);hold on;
xlabel('y_1')
ylabel('y_2')
lags = [1]; % this is where we specify the vector of taus

sol = dde23(@ddefun,lags,@history,tspan)

y1 = sol.y(1,:); % note the y-solution is a row-wise matrix.
y2 = sol.y(2,:);
plot(y1,y2);

function dYdt = myode(t,Y)
% ordinary model
y1 = Y(1);
y2 = Y(2);

dy1dt = a*y1 + b*y1*y2;
dy2dt = c*y2 + d*y1*y2;

dYdt = [dy1dt; dy2dt];
end


% we need a history function that is a column vector giving the value of y1
% and y2 in the past. Here we make these constant values.
function y = history(t)
y = [80; 30];
end

% we define the function for the delay. the Y variable is the same as you
% should be used to from an ordinary differential equation. Z is the values
% of all the delayed variables.
function  dYdt = ddefun(t,Y,Z)
m = 200; % additional variable
y1 = Y(1);
y2 = Y(2);

% Z(:,1) = [y1(t - tau_1); y2(t - tau_2)]
y1_tau1 = Z(1,1);
y2_tau1 = Z(2,1);

dy1dt = a*y1*(1-y1/m) + b*y1*y2;
dy2dt = c*y2 + d*y1_tau1*y2_tau1;

dYdt = [dy1dt; dy2dt];
end



end