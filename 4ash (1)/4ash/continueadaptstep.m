






for l = 1:M.np

        M.SS(l).SS(:,1) = M.SS(l).SS(:,k-1);
        M.GV(l,:,1)     = M.GV(l,:,k-1);

end

%%

k=1;

while k < M.opt.N %&& r(k)>10*tol % integrate through time
        
        % --- Prediction step ---
%         dt = mindt;
        k=k+1;
        for l = 1:M.nc % Compute recieved conductances
            Gout(l,:) = P.A(l).M*squeeze(M.GV(:,l,k-1)) + P.C(l).M;
        end
        
        for l = 1:M.np   % cycle through populations
            
            P1      = M.P(l);
            SSS2    = [reshape(M.SS(l).SS(:,k-1),[],1,1); squeeze(M.GV(l,:,k-1))'];
            [Qdpdt] = fx_LIFpopME(Gout(:,l),P1);
            M1      = expm(.5*dt*Qdpdt);
            SSS     = M1*SSS2; % Midpoint
            SSS1    = M1*SSS;  % Endpoint
            
            M.SS(l).SS(:,k) = SSS(1:P1.LVV);
            M.SS0(l).SS = SSS1(1:P1.LVV);
            
            M.GV(l,:,k)     = SSS(P1.LVV+1:end);
            M.GV0(l,:)      = SSS1(P1.LVV+1:end);
            
        end
        
        % --- Correction step ---9
        
        
        for l = 1:M.nc % Compute recieved conductances
            Gout(l,:) = P.A(l).M*(squeeze(M.GV(:,l,k))) + P.C(l).M; % Midpoint
        end
        
        for l = 1:M.np   % cycle through populationsexp(-6)
            P1      = M.P(l);
            SSS2    = [reshape(M.SS(l).SS(:,k-1),[],1,1); squeeze(M.GV(l,:,k-1))'];
            [Qdpdt] = fx_LIFpopME(Gout(:,l),P1);
            SSS     = expm(dt*Qdpdt)*SSS2;
            
            M.SS(l).SS(:,k) = SSS(1:P1.LVV);
            M.GV(l,:,k)     = SSS(P1.LVV+1:end);
        end
        
        M.opt.t(k) = M.opt.t(k-1) + dt;
        
        D = (M.GV0 - squeeze(M.GV(:,:,k))).^2;
        D = sum(D(:));
        
        for l=1:M.np
            D = D + sum((M.SS0(l).SS - squeeze(M.SS(l).SS(:,k))).^2);
        end
        
        
%         D = 0;
        if D >tol
            dt = dt/2;
            k  = k - 1;
%             disp(['timestep -: ' num2str(dt)])
        elseif (sum(D(:)) < tol/8) %&& dt<M.opt.dt*2
            dt = 2*dt;
            inlineprogress(k,M.opt.N)
%             disp(['timestep +: ' num2str(dt)])
        else
            inlineprogress(k,M.opt.N)
        end
        
        if ~(sum(D(:)) > tol)
            
            if k > 10
                [r0, t0] = convergence_check(M.GV,M.opt.t,k,M.opt.T0);
                r(k) = r0;
            end
            % -----------------------
            % --- LFP computation ---
            % -----------------------
            
            if ~isfield(M.opt,'wLFP'); M.opt.wLFP = ones(size(M.P(1).T)); end;
            if ~isfield(M.opt,'LocalElectrode'); M.opt.LocalElectrode = 0; end;
            
            for l = 1:M.np
                P1 = M.P(l);
                FV = P1.FvarV;
                FV(P1.Tr+1:end,:) = FV(repmat(P1.R,P1.LQ,1),:); % heuristic on refractory period
                M.LFP(l,k) = M.opt.wLFP*(Gout(:,l).*(FV'*squeeze(M.SS(l).SS(:,k)*P1.Vres))); % [nS*mV] = [pA] in pico Ampere (per neuron)
                if M.opt.LocalElectrode
                    M.LFP(l,k) = M.LFP(l,k) + squeeze(M.SS(l).SS(P1.Tr+1,k)*P1.Vres)*P1.C*(P1.Vt-P1.Vr)/.001; % Current for repolarization in 1 ms (per neuron)
                end
            end
            
            % --- Current computation ---
            for m= 1:M.np
                P1  =  M.P(m);
                for l = 1:M.nc
                    M.Currents(m,l,k) = Gout(l,m)*(P1.Vg(1)-P1.Vg(2)) + P1.gl*(P1.Vl-P1.Vg(2)); % [nS*mV] = [pA] in pico Ampere (per neuron)
                end
            end
            
            
            if ~isfield(M.opt,'pop'); M.opt.pop = 2; end;
            pop = M.opt.pop;
            
            if M.opt.Drawplots
                
                for l = 1:M.np   % cycle through populationsexp(-6)
                    P1 = M.P(l);
                    subplot(M.np+2,1,l)
                    plot(P1.VV,(M.SS(l).SS(:,k)));axis([P1.VV(1) P1.VV(end) -.001 .2]);
                    drawnow;
                end
                
                
                subplot(M.np+2,1,M.np+1)
                plot(M.opt.t,M.LFP(pop,:)','k');
                hold on
                plot(M.opt.t,squeeze(M.Currents(pop,:,:))','.','MarkerSize',1)
                hold off
                 
                subplot(M.np+2,1,M.np+2)
%                 plot(squeeze(M.Currents(pop,1,:)),squeeze(M.Currents(pop,2,:)),'.','MarkerSize',1);
                plot(M.opt.t,r,'.'); hold on;plot(t0,r0,'r.'); hold off;
                xlabel('E (pA)'); %xlim([for l = 1:M.np
                ylabel('I (pA)'); %ylim([-100 100]+EIcurrents(k,2));
                drawnow;
            end
        end
        
        
end
    
    


%%
M2 = M; mindt = min(diff(M.opt.t));

%%
figure
% plot(M2.opt.t,squeeze(M2.Currents(pop,:,:))');%,'-','MarkerSize',1)
hold on
plot(M.opt.t,squeeze(M.Currents(pop,:,:))','.','MarkerSize',1)
%%

figure

derivative = bsxfun(@ldivide,diff(M.opt.t)',diff(squeeze(M.Currents(pop,:,:)),1,2));
plot((M.opt.t(1:end-1)),derivative,'-.','MarkerSize',1)




