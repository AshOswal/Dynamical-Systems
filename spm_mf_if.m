function spm_mf_if

% model dynamics
%==========================================================================
 
% create inputs
%--------------------------------------------------------------------------
dt      = 1e-4;
n       = 100e-3/dt;
u       = 8;
t       = (1:n)*dt*1e3;
input   = simulate_poisson_train(600,600,1e-3,1,1e-3);
x0      = [-80; 0.01];

% integrate with input - output = E{x} and create response
%--------------------------------------------------------------------------
for j = 1:32
    x = x0;
    for i = 1:n
        dx      = spm_fx_if([sum(input(:,n)); x]) + 1*randn(size(x0));
        x       = x + dx*dt;
        V(i,j)  = x(1); % state space voltage
        T(i,j)  = x(2); % state space time
    end
end

% single trajectories over time
%--------------------------------------------------------------------------
%spm_figure('GetWin','Figure 1'); clf
figure;
subplot(2,2,1)
plot(t,V)
title('Trajectories over time','FontSize',16)
xlabel('time (ms)','FontSize',12)
ylabel('depolarization','FontSize',12)
axis square
grid on
 
clear mov
subplot(2,2,2)
V = V(1:2:end,:);
T = T(1:2:end,:);
for i = 1:length(V)
    
    % plot
    %----------------------------------------------------------------------
    plot(T(i,:)*1000,V(i,:),'o','MarkerSize',8)
    axis([0 32 -100 -40])
    title({'Trajectories in state-space';'Click axis for movie'},'FontSize',16)
    xlabel('peri-spike time (ms)','FontSize',12)
    ylabel('depolarization','FontSize',12)
    axis square
    grid on
    drawnow
    
    % addframe
    %----------------------------------------------------------------------
    mov(i) = getframe(gca);
    
end
 
% set ButtonDownFcn
%--------------------------------------------------------------------------
set(gca,'Userdata',{mov,16})
%set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
 
 




