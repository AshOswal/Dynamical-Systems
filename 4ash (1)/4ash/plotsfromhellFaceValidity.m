%% plotsfromhellFaceValidity


testrun=[2 2 3 2 2 2 3 3 3 3];
for n=7%1:10;
    
    figure(n)
    
    load(['FaceFR' num2str(testrun(n)) '/M' num2str(n) '.mat'])
    load(['randYfr/Y' num2str(n) '.mat'])
    
    
    
    if isfield(M,'f')
        f = M.f;
    else
        [dfdp,f] = spm_DCM_lifpopsys_LC_adapstep_y_par_fr(M.Ep,M,0*spm_vec(M.pE));
        f= spm_unvec(spm_vec(f)*0,f);
    end
   
    dp = 0;
    k=2;
    Cm = lines(7); Cm(3,:) = [];
    fp = f;
    y = Y.y;
    imind=[1 3 2 4];
    
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
    
end
%%

p = zeros(10,1);

for n=1:10;
    
    load(['FaceFR' num2str(testrun(n)) '/M' num2str(n) '.mat'])
    load(['randYfr/Y' num2str(n) '.mat'])
    
    figure(n)
%     plot(spm_vec(M.pE),'*')
%     hold on
%     plot(spm_vec(M.pE) + 2*sqrt(spm_vec(M.pC)),'.')
%     plot(spm_vec(M.pE) - 2*sqrt(spm_vec(M.pC)),'.')
%     
%     plot(spm_vec(M.Ep),'k*')
%     plot(spm_vec(M.Ep) + 2*sqrt(diag(M.Cp)),'k.')
%     plot(spm_vec(M.Ep) - 2*sqrt(diag(M.Cp)),'k.')
%     
%     plot(spm_vec(P),'r*')
    
    prior = spm_vec(M.pE);
    sigprior = spm_vec(M.pC);
    mu = spm_vec(M.Ep);
    x  = spm_vec(P);
    sig = M.Cp;
    ind = abs(diag(sig))>10^-3;
    prior = prior(ind);
    sigprior = sigprior(ind);
    mu = mu(ind);
    x  = x(ind);
    sig = sig(ind,ind);
    p(n) = 1-chi2cdf((mu-x)'/sig*(mu-x),length(x));
    
    plot(0*prior,'*')
    hold on
    plot( 2*ones(size(sigprior)),'.')
    plot(-2*ones(size(sigprior)),'.')
    
    plot((mu-prior)./sqrt(sigprior),'k*')
    plot((mu-prior + 2*sqrt(diag(sig)))./sqrt(sigprior),'k.')
    plot((mu-prior - 2*sqrt(diag(sig)))./sqrt(sigprior),'k.')
    
    plot((x-prior)./sqrt(sigprior),'r*')
    
grid on
    
end
