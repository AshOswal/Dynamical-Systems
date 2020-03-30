x = [-10:.01:10]';
Daw  = mfun('dawson',x);
Hip  = 1./(2*x);
F    = sqrt(pi)/2*sqrt(1 - exp((4/pi)*x.^2)).*sign(x)/1i;
S    = x - 2/3*x.^3 + 4/15*x.^5 - 8/105*x.^7 + 16/945*x.^9; % series arround 0
S2   = 1./(2*x) + 1./(4*x.^3); + 3./(8*x.^5) + 15./(16*x.^7) + 105./(32*x.^9);
L1   = 1./(2+4*x.^2);
L2   = 1./(12*x + 8*x.^3);
L4   = 1./(12+48*x.^2+16*x.^4);
dF   = 2*x.*Daw + exp(-2*x.^2);

plot(x,[Daw exp(-x.^2).*F S2 ])

axis([-10 10 -1 1])

%%
N     = 20;
coef  = zeros(N,1);
coefi = zeros(N,1);
X     = zeros(length(x),N);
Xi    = zeros(length(x),N);
n     =  (1:2:2*N+1)';
% Daw = mfun('dawson',x);
for k = 1:N
    coef(k)  = (-1)^(k-1)*2^(k-1)/prod(n(1:k));
    coefi(k) = prod([1; n(1:k-1)])/2^(k);
    X(:,k)   = x.^n(k); 
    Xi(:,k)  = (x).^-(n(k)); 
end


A = X*coef;
coefM = repmat(coefi,1,N);
coefM = triu(coefM);
B = Xi*coefM;
plot(x,[Daw A B])
% plot(x,[(repmat(A,1,N)-B),0*x]);plot(x,[A,B,0*x])
axis([-10 10 -1 1])

%%
P02 = -Daw(1:end/2);
P01 = exp(-x(1:end/6).^2);
P01 = P01/P01(end)*P02(1);
P03 = P02(1)*2*ones(length(x)/6,1);
P0 = [P01;P02;P03];
P0 = P0/norm(P0);
plot(P0)

x1 = -length(P0)/2:length(P0)/2-1;
x1 = x1'/norm(x1);
x1 = x1.*P0;
P1 = x1-(x1'*P0)*P0;
P1 = P1/norm(P1);
plot([P0 P1])


