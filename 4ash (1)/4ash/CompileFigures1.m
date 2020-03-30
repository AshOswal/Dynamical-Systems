%% Compile Figures


figure('color','w')

load('LIF_FP_500mus_250muV_res_steady.mat')
SSoverT2 = SSoverT;
SSoverT(:,1,:)=SSoverT(:,2,:);
SSoverT(:,2,:)=SSoverT2(:,1,:);
mS = max(max(SSoverT(P.Tr+1,:,:)*P.Fsteady(P.Tr+1)));
mS = 1000;
k = N;
 for l = 1:np
     subplot(np+1,3,3*(l-1)+1);plot(P.VV,(SSoverT(:,l,k)));axis([P.Vmin P.Vmax -.001 .2]);
 end
 subplot(np+1,3,np*3+1);
 plot(dt:dt:N*dt,squeeze(SSoverT(P.Tr+1,:,:)*P.Fsteady(P.Tr+1))')
 ylim([0 mS])
 xlim([dt N*dt])
 hold on
 plot([dt*k dt*k],[0 mS],'k','LineWidth',2)
 hold off
 drawnow;
 
load('LIF_FP_500mus_250muV_res_spikes_6hz.mat')
SSoverT2 = SSoverT;
SSoverT(:,1,:)=SSoverT(:,2,:);
SSoverT(:,2,:)=SSoverT2(:,1,:);
mS = max(max(SSoverT(P.Tr+1,:,:)*P.Fsteady(P.Tr+1)));
mS = 1000;
k = N;
 for l = 1:np
     subplot(np+1,3,3*(l-1)+2);plot(P.VV,(SSoverT(:,l,k)));axis([P.Vmin P.Vmax -.001 .2]);
 end
 subplot(np+1,3,np*3+2);
 plot(dt:dt:N*dt,squeeze(SSoverT(P.Tr+1,:,:)*P.Fsteady(P.Tr+1))')
 ylim([0 mS])
 xlim([dt N*dt])
 hold on
 plot([dt*k dt*k],[0 mS],'k','LineWidth',2)
 hold off
 drawnow;
 
load('LIF_FP_500mus_250muV_res_spikes.mat')
SSoverT2 = SSoverT;
SSoverT(:,1,:)=SSoverT(:,2,:);
SSoverT(:,2,:)=SSoverT2(:,1,:);

mS = max(max(SSoverT(P.Tr+1,:,:)*P.Fsteady(P.Tr+1)));
mS = 1000;
k = N;
 for l = 1:np
     subplot(np+1,3,3*(l-1)+3);plot(P.VV,(SSoverT(:,l,k)));axis([P.Vmin P.Vmax -.001 .2]);
 end
 subplot(np+1,3,np*3+3);
 plot(dt:dt:N*dt,squeeze(SSoverT(P.Tr+1,:,:)*P.Fsteady(P.Tr+1))')
 ylim([0 mS])
 xlim([dt N*dt])
 hold on
 plot([dt*k dt*k],[0 mS],'k','LineWidth',2)
 hold off
 drawnow;

%%
subplot(np+1,3,1)
subplot(np+1,3,2)
title('\fontsize{14} Inhibitory Interneurons','Color','r')
xlabel('State Space (mV)')
subplot(np+1,3,3)
subplot(np+1,3,4)
ylabel('\fontsize{14} Probability density')
subplot(np+1,3,5)
title('\fontsize{14} Excitatory Spiny Cells','Color',[0 .5 0])
xlabel('State Space (mV)')
subplot(np+1,3,6)
subplot(np+1,3,7)
subplot(np+1,3,8)
title('\fontsize{14} Excitatory Pyramidal Cells','Color','b')
xlabel('State Space (mV)')
subplot(np+1,3,9)
subplot(np+1,3,10)
ylabel({'\fontsize{12} Firing Rate';'(spikes s^-^1 neuron^-^1)'})
subplot(np+1,3,11)
xlabel('Time (s)')
subplot(np+1,3,12)

