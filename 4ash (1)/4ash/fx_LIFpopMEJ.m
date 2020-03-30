function [f,J] = fx_LIFpopMEJ(NodeP,G,U,S)
% Full flow for a system of interconnected LIFpop
% G(:,:,1) = GE; G(:,:,2) = GI; U =  U(:,1)*C'; S = [SSoverT(:,:,kini-1); squeeze(GV(:,:,kini-1))];
[ns np] = size(S);
ng      = size(G,3);
f       = zeros(np*ns,1);
J       = zeros(np*ns);
Go      = zeros(ng,np);

for g = 1:ng
    Go(g,:) = G(:,:,g)*S(ns-ng+g,:)' + U(g,:)';
end

for k = 1:np
    J((k-1)*ns+1:k*ns,(k-1)*ns+1:k*ns) = ...
        J((k-1)*ns+1:k*ns,(k-1)*ns+1:k*ns) + fx_LIFpopME(Go(:,k),NodeP(k));
end

f = J*S(:);
for k = 1:np    
    for l = 1:np
        
        for g = 1:ng
            dJ0  = (G(l,k,g)/NodeP(l).Vres)*NodeP(l).FvarV(:,g).*S((l-1)*ns+1:(l*ns-ng))';
            J((l-1)*ns+1:(l*ns-ng),(k*ns-ng+g)) = J((l-1)*ns+1:(l*ns-ng),(k*ns-ng+g))...
                + [dJ0(2:end); 0] - [-dJ0(1); dJ0(1:end-1)];
        end
    end
end


return


%% Check eigen values
[V S1] = eig(full(J));
[B,IX] = sort(real(diag(S1)));

figure
subplot(2,1,1)
plot(real(diag(S1)),imag(diag(S1)),'.')
plot(B)
% xlim([-200,50,])
subplot(2,1,2)
plot(real(V(:,IX(end-1:end))))

%%
S2 = S;
J2       = zeros(np*ns);

for g = 1:ng
    Go2(g,:) = G(:,:,g)*S2(ns-ng+g,:)' + U(g,:)';
end

for k = 1:np
    J2((k-1)*ns+1:k*ns,(k-1)*ns+1:k*ns) = ...
        J2((k-1)*ns+1:k*ns,(k-1)*ns+1:k*ns) + fx_LIFpopME(Go2(:,k),NodeP(k));
end

f1 = J2*S2(:);

deltaE = +.0001;
stateE = 102;
S2(stateE)=S2(stateE) + deltaE;
J2       = zeros(np*ns);

for g = 1:ng
    Go2(g,:) = G(:,:,g)*S2(ns-ng+g,:)' + U(g,:)';
end

for k = 1:np
    J2((k-1)*ns+1:k*ns,(k-1)*ns+1:k*ns) = ...
        J2((k-1)*ns+1:k*ns,(k-1)*ns+1:k*ns) + fx_LIFpopME(Go2(:,k),NodeP(k));
end

f2 = J2*S2(:);
plot([(f2-f1)/deltaE J(:,stateE)])
end
