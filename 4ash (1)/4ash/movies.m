%% Film capture

SSoverT2 = SSoverT;
SSoverT(:,1,:)=SSoverT(:,2,:);
SSoverT(:,2,:)=SSoverT2(:,1,:);

figure('color','w')
mS = max(max(SSoverT(P.Tr+1,:,:)*P.Fsteady(P.Tr+1)));
for k = 1:N
 for l = 1:np
     subplot(np+1,1,l);plot(P.VV,(SSoverT(:,l,k)));axis([P.Vmin P.Vres*P.LVV+P.Vmin -.001 .2]);
 end
 subplot(np+1,1,np+1);
 plot(dt:dt:N*dt,squeeze(SSoverT(P.Tr+1,:,:)*P.Fsteady(P.Tr+1))')
 ylim([0 mS])
 hold on
 plot([dt*k dt*k],[0 mS],'k','LineWidth',2)
 hold off
 drawnow;
 F(k) = getframe(gcf);
end


%%
writerObj = VideoWriter('test2.mp4','MPEG-4');
open(writerObj);
for k = 1:N
    writeVideo(writerObj,F(k));
    disp(k)
end
close(writerObj);


%%

figure('color','w')

%%
mS = max(max(SSoverT(P.Tr+1,:,:)*P.Fsteady(P.Tr+1)));
for k = 1:N
 for l = 1:np
     subplot(np+1,1,l);plot(P.VV,(SSoverT(:,l,k)));axis([P.Vmin P.Vres*P.LVV+P.Vmin -.001 .2]);
     if l == pp; hold on; plot(VL(10*k,1:6),.01*ones(1,6),'ro'); hold off; end
 end
 subplot(np+1,1,np+1);
 plot(dt:dt:10*k*dt,VL(1:10*k,1:6))
 ylim([-90 30])
 xlim([dt 10000*dt])
 hold on
 plot([10*dt*k 10*dt*k],[-90 0],'k','LineWidth',2)
 plot(dt:10*dt:10*k*dt,squeeze(.05*LFP(2,1:k)))
 hold off
 
subplot(np+1,1,2)
ylabel('Probability density')
subplot(np+1,1,np)
xlabel('State Space ("mV")')
subplot(np+1,1,np+1)
xlabel('Time (s)')
ylabel({'Firing Rate';'(spikes s^-^1)'})

subplot(np+1,1,1)
title('Inhibitory Interneurons')
subplot(np+1,1,2)
title('Excitatory Pyramidal Cells')
subplot(np+1,1,3)
title('Local Field Potential')


drawnow;
F(k) = getframe(gcf);

 
end


%%

subplot(np+1,1,2)
ylabel('Probability density')
subplot(np+1,1,np)
xlabel('State Space ("mV")')
subplot(np+1,1,np+1)
xlabel('Time (s)')
ylabel({'Firing Rate';'(spikes s^-^1)'})

subplot(np+1,1,1)
title('Inhibitory Interneurons')
subplot(np+1,1,2)
title('Excitatory Spiny Cells')
subplot(np+1,1,3)
title('Excitatory Pyramidal Cells')

