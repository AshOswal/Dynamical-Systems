%% Compile Figures


figure('color','w')

load('LIF_FP_500mus_250muV_res_spikes_sparse2.mat')
mS = max(max(SSoverT(P.Tr+1,:,:)*P.Fsteady(P.Tr+1)));
mS = 200;
k = N;
 for l = 1:np
     subplot(np+1,1,l);plot(P.VV,(SSoverT(:,l,k)));axis([P.Vmin P.VV(end) -.001 .2]);
 end
 subplot(np+1,1,np+1);
 plot(dt:dt:N*dt,squeeze(SSoverT(P.Tr+1,:,:)*P.Fsteady(P.Tr+1))')
 ylim([0 mS])
 xlim([dt N*dt])
 hold on
 plot([dt*k dt*k],[0 mS],'k','LineWidth',2)
 hold off
 drawnow;


%%
subplot(np+1,1,1)
title('\fontsize{14} Inhibitory Interneurons','Color','b')
xlabel('State Space (mV)')
subplot(np+1,1,2)
title('\fontsize{14} Excitatory Pyramidal Cells','Color',[0 .5 0])
xlabel('State Space (mV)')
ylabel('\fontsize{14} Probability density')
subplot(np+1,1,3)
ylabel({'\fontsize{12} Firing Rate';'(spikes s^-^1 neuron^-^1)'})
xlabel('Time (s)')

