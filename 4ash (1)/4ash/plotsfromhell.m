%% plotsfromhell


if isfield(M,'f')
    f = M.f;
else
    [dfdp,f] = spm_DCM_lifpopsys_LC_adapstep_y_par_fr(M.Ep,M,0*spm_vec(M.pE));
    f= spm_unvec(spm_vec(f)*0,f);
end
%%

%   figure
%   f = x;
%   dfdp=dxdp;
  dp = 0;
  k=2;
    Cm = lines(7); Cm(3,:) = [];
    fp = f;
    y = Y.y;
       imind=[1 3 2 4];
%     
%     for k1 = 1:4
%         subplot(3,2,imind(k1))
%         hold on
%         plot(y(:,k1),'r','LineWidth',2)
%         plot(f.lc(:,k1),'Color',Cm(mod(k,6)+1,:))
%         plot(fp.lc(:,k1),'.','Color',Cm(mod(k,6)+1,:),'MarkerSize',1)
%         hold off
%     end
%     
%     subplot(3,2,5)
%     hold on
%     plot(-y(:,2),y(:,1),'r','LineWidth',2)
%     plot(-f.lc(:,2),f.lc(:,1),'Color',Cm(mod(k,6)+1,:))
%     plot(-fp.lc(:,2),fp.lc(:,1),'.','Color',Cm(mod(k,6)+1,:),'MarkerSize',1)
%     title(['Efr = ' num2str(f.fr(1,1)) '  dfr = ' num2str(fp.fr(1,1)) '  fr = ' num2str(y(1,end-2))])
%     xlabel(['Ec = ' num2str(f.c(1)) '  df = ' num2str(fp.c(1)) '  f = ' num2str(y(1,end))])
%      
%     subplot(3,2,6)
%     hold on
%     plot(-y(:,4),y(:,3),'r','LineWidth',2)
%     plot(-f.lc(:,4),f.lc(:,3),'Color',Cm(mod(k,6)+1,:))
%     plot(-fp.lc(:,4),fp.lc(:,3),'.','Color',Cm(mod(k,6)+1,:),'MarkerSize',1)
% %     hold off
%     title(['Efr = ' num2str(f.fr(1,2)) '  dfr = ' num2str(fp.fr(1,2)) '  fr = ' num2str(y(1,end-1))])
%     xlabel(['Ef = ' num2str(f.f(1)) '  df = ' num2str(fp.f(1)) '  f = ' num2str(y(1,end-3))])
%     drawnow
    
    figure
    
    for k1 = 1:4
        subplot(4,2,imind(k1))
        hold on
        plot(y(:,k1),'r','LineWidth',2)
        plot(f.lc(:,k1),'Color',Cm(mod(k,6)+1,:))
        plot(fp.lc(:,k1),'.','Color',Cm(mod(k,6)+1,:),'MarkerSize',1)
        hold off
    end
    
    subplot(4,2,5)
    hold on
    plot(y(:,6),'r','LineWidth',2)
    plot(f.fr(:,1),'Color',Cm(mod(k,6)+1,:))
    plot(fp.fr(:,1),'.','Color',Cm(mod(k,6)+1,:),'MarkerSize',1)
    hold off
    
    subplot(4,2,6)
    hold on
    plot(y(:,7),'r','LineWidth',2)
    plot(f.fr(:,2),'Color',Cm(mod(k,6)+1,:))
    plot(fp.fr(:,2),'.','Color',Cm(mod(k,6)+1,:),'MarkerSize',1)
    hold off
    
    subplot(4,2,7)
    hold on
    plot(-y(:,2),y(:,1),'r','LineWidth',2)
    plot(-f.lc(:,2),f.lc(:,1),'Color',Cm(mod(k,6)+1,:))
    plot(-fp.lc(:,2),fp.lc(:,1),'.','Color',Cm(mod(k,6)+1,:),'MarkerSize',1)
    title(['Efr = ' num2str(f.fr(1,1)) '  dfr = ' num2str(fp.fr(1,1)) '  fr = ' num2str(y(1,end-2))])
    xlabel(['Ec = ' num2str(f.c(1)) '  dc = ' num2str(fp.c(1)) '  c = ' num2str(y(1,end))])
     
    subplot(4,2,8)
    hold on
    plot(-y(:,4),y(:,3),'r','LineWidth',2)
    plot(-f.lc(:,4),f.lc(:,3),'Color',Cm(mod(k,6)+1,:))
    plot(-fp.lc(:,4),fp.lc(:,3),'.','Color',Cm(mod(k,6)+1,:),'MarkerSize',1)
%     hold off
    title(['Efr = ' num2str(f.fr(1,2)) '  dfr = ' num2str(fp.fr(1,2)) '  fr = ' num2str(y(1,end-1))])
    xlabel(['Ef = ' num2str(f.f(1)) '  df = ' num2str(fp.f(1)) '  f = ' num2str(y(1,end-3))])
    drawnow