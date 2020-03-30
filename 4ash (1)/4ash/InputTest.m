% load Ep
TestParam8CA1;
load('Munt4_1')

Ep0 = M.Ep;
Ep = Ep0;
Ninput = 40;
Flist = zeros(Ninput,1);
Alist = zeros(Ninput,1);
flist = cell(Ninput,1);
Drawplots = 0;
for inputU = 1:Ninput
    fprintf('Input %i', inputU)
%     Ep.A(1).M(1,2) = -32; 
%     Ep.A(3).M(2,2) = -32; 
    
    Ep.C(1).M   = Ep0.C(1).M + Ninput*.005*(inputU-Ninput/2); 
    
    [dfdp,f] = spm_DCM_lifpopsys_LC_adapstep_y_par_fr(Ep,M,0*spm_vec(M.pE));
%     if ~isempty(M.f)
%         f = M.f;
%     else
%         f= spm_unvec(spm_vec(f)*0,f);
%     end
    
    dp = 0;
    k=2;
    Cm = lines(7); Cm(3,:) = [];
    fp = f;
    y = Y.y;
    Flist(inputU) = f.f(1);
    Alist(inputU) = norm(f.lc(:,1));
    imind=[1 3 2 4];
    if Drawplots
        figure(inputU)
        
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
        %     title(['Efr = ' num2str(f.fr(1,1)) '  dfr = ' num2str(fp.fr(1,1)) '  fr = ' num2str(y(1,end-2))])
        xlabel(['Ec = ' num2str(f.c(1)) '  dc = ' num2str(fp.c(1)) '  c = ' num2str(y(1,end))])
        
        subplot(4,2,8)
        hold on
        plot(-y(:,4),y(:,3),'r','LineWidth',2)
        plot(-f.lc(:,4),f.lc(:,3),'Color',Cm(mod(k,6)+1,:))
        plot(-fp.lc(:,4),fp.lc(:,3),'.','Color',Cm(mod(k,6)+1,:),'MarkerSize',1)
        %     hold off
        %     title(['Efr = ' num2str(f.fr(1,2)) '  dfr = ' num2str(fp.fr(1,2)) '  fr = ' num2str(y(1,end-1))])
        xlabel(['Ef = ' num2str(f.f(1)) '  df = ' num2str(fp.f(1)) '  f = ' num2str(y(1,end-3))])
        drawnow
    end
    flist{inputU} = f;
end




%%

for inputU = 1:Ninput
    
    f = flist{inputU};
    Flist(inputU) = f.f(1);
    Alist(inputU) = max(f.lc(:,1)) - min(f.lc(:,1)) ;
    
end
    





