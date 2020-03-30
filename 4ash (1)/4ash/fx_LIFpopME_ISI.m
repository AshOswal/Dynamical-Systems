function [Qdpdt] = fx_LIFpopME_ISI(Go,P)

% Some code optimization is still possible, profiting by saving
% Diff matrices from iteration to iteration.

m = P.LVV;
n = length(P.G);
T = P.Tr;
R = P.R;
Fsteady = P.Fsteady;
FvarV   = P.FvarV;

FvarG = zeros(size(Fsteady)); 
for k = 1:n
    FvarG(:) = FvarG(:) - (Go(k)*P.G(k) + P.DGL(k)*P.DL(k))*FvarV(:,k);
end

F = (Fsteady + FvarG);

FS1  =  sparse(1:m,1:m,F(:),m,m);
D0   =  spdiags(ones(m,1),0,m,m);
D1   =  spdiags(ones(m,1),1,m,m);
D_1  =  spdiags(ones(m,1),-1,m,m);

Flow          = (-D0*abs(FS1) - D1*(FS1.*(FS1<0)) + D_1*(FS1.*(FS1>0)))/P.Vres;
Flow(1,1)     = -Flow(2,1);

D2              = spdiags([ones(m,1) -2*ones(m,1) ones(m,1)],[-1 0 1],m,m);%LIFDiff21M(m,1);
D2(T,T-2:T+2)   = 2*[0 1/2 -1 0 0];   %2*[1/2 -1 1/2 0 0];
D2(T-1,T-3:T+1) = 2*[0 1/2 -1 1/2 0]; %2*[0 1/2 -1 1/2 0];
D2(1,1)         = -1;
% D2(T-2,T-4:T)   = 2*[-1/24 2/3 -5/4 2/3 0];
D2(T+1:end,:)   = 0;
D2(T+1,T)       = 1;

% noise for all channel families

DS1 = zeros(m,m);

for k = 1:length(Go)
    DS1 = DS1 + sparse(1:m,1:m, P.DGL(k)*.5*(P.DL(k)^2)*abs(FvarV(:,k)) ,m,m);
end


Diffusion = 1/2*D2*(P.D*speye(m) + DS1)/(P.Vres^2);

% Diffusion(T,T) = 0;

ResetOperator = FS1(end,:)/P.Vres;

Qdpdt         = Flow + Diffusion;
% Qdpdt(R,:)    = Qdpdt(R,:)   + ResetOperator;

% SrateOperator and Gi flow
Qdpdt(m+1:m+n,P.delay) = Qdpdt(P.delay+1,P.delay).*(1./P.T)*P.Vres;
Qdpdt(m+1:m+n,m+1:m+n) = -diag((1./P.T));

% Population fiering rate

% dGdt  = (1./P.T).*((1-Gi(:)').*(P.p+Srate)-Gi(:)'); % MLGFlow(G,Srate,P);

end

% %
% subplot(2,1,1)
% imagesc(Qdpdt)
% subplot(2,1,2)
% A2 = null(full(Qdpdt));
% A = A2; A(T)=0;
% plot([Qdpdt*A A A2])
% % 
% 
% [V S] = eig(full(Qdpdt));
% [B,IX] = sort(real(diag(S)));
% subplot(2,1,1)
% plot(B)
% subplot(2,1,2)
% plot(real([V(:,IX(end))]))
% 
% S1 = S;
% S1(IX(end),IX(end))=0;
% imagesc(real(V*S1*inv(V)))
% imagesc(real(V*S1*inv(V))-Qdpdt)

% % 
% % % 
% plot(real([E(:,1) ]))
