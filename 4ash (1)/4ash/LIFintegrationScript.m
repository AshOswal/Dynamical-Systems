%% SimLIFneuron
% 
VL = []; FR = zeros(1,10*N);
N2 = 10;
VL = [VL zeros(10*N,N2)];

pp = 2;
dt1 = dt/10;

for k2 = 1:N2
disp(k2)
    
A = SSoverT(:,l,1);
B = cumsum(A);
V = interp1(B(B<.9999 & B>.0001),P.VV(B<.9999 & B>.0001),rand);
VV = zeros(1,10*N);

if V > P.Vt
    k = round(P.Rt/dt1*(V-P.Vt)/(P.VV(end)-P.Vt))+2;
    V = P.Vr;
    VV(k-1) = V;
else
    VV(1) = V;
    k = 2;
end

while k <= 10*N

    Gout(1,:) = GE*GV(1,:,ceil((k)/10))'+ C*U(ceil((k)/10));
    Gout(2,:) = GI*GV(2,:,ceil((k)/10))';
%   Gout = .4*ones(size(Gout));

    r = sqrt(P.D*dt1)*randn;
    dVdt  = fx_LIF(VV(k-1),Gout(:,pp),0,0,P);
    dVdt2 = fx_LIF(VV(k-1) + dt1*dVdt + r ,Gout(:,pp),0,0,P);
    VV(k) = VV(k-1) + .5*dt1*(dVdt + dVdt2) + r;
    if VV(k)>P.Vt
        dtFire = (VV(k)-P.Vt)/(VV(k)-VV(k-1));
        VV(k) = 0;
        FR(k) = FR(k)+1;
        k = k + round(P.Rt/dt1) - 1;
        if k <10*N
            r = sqrt(P.D*dtFire*dt1)*randn;
            dVdt  = fx_LIF(P.Vr,Gout(:,pp),0,0,P);
            dVdt2 = fx_LIF(P.Vr + dt1*dtFire*dVdt + r ,Gout(:,pp),0,0,P);
            VV(k) = P.Vr + .5*dt1*dtFire*(dVdt + dVdt2) + r;
        end
    end
    k = k+1;
end

% plot(VV)

VL(:,N2+1-k2) = VV';

end

% FR = sum(VL==0,2);



%%
ttt = dt1:dt1:10*N*dt1;
Kernel = exp(-1/P.T(1)*ttt);
Kernel = Kernel/sum(Kernel);
FRC = conv(FR,Kernel);
FRC = FRC/size(VL,2)/dt1;
FRCD = FRC(10:10:10*N);
figure
plot([FRCD' squeeze(GV(1,pp,:))])

%%

% hist(VL(end-1,:),P.VV)
j = 478;
AAA = VL(end-10*(j+1):end-10*j,:);
HHH = hist(AAA(:),P.VV);
HHH = HHH/sum(HHH);
plot(P.VV,[1/P.Vres*HHH' SSoverT(:,pp,end-j)])
axis([-100 -30 0 .2])

%% PSTH

[row,col] = find(VL(end/4:end-1,:)==0);
Drow=diff(row);
Drow(Drow<=1) = [];
Drow = Drow*dt1;
% figure
ISIh = hist(Drow,(1:300)/2000);
ISIh = ISIh/sum(ISIh);

%% stationary solution
StSol = LIF_StSol(Gout(:,pp),P,0);
StSol = StSol/sum(StSol);
j = 0;
AAA = VL(end-10*(j+150):end-10*j,:);
HHH = hist(AAA(:),P.VV);
HHH = HHH/sum(HHH);
plot(P.VV,[HHH' StSol'])
axis([-100 -30 0 .02])

%%
Nfr = 500;
StFR0 = zeros(P.LVV,Nfr);
StFR1 = zeros(P.LVV,Nfr);
StFR4pi= zeros(P.LVV,Nfr);
StFR7= zeros(P.LVV,Nfr);
for k =1:Nfr
    StFR0(:,k) = LIF_StSol([2*k/Nfr;0],P,0);
%     StFR1(:,k) = LIF_StSol([2*k/Nfr;0],P,1);
%     StFR4pi(:,k) = LIF_StSol([2*k/Nfr;0],P,4/pi);
%     StFR7(:,k) = LIF_StSol([2*k/Nfr;0],P,.736);
%     StFR8(:,k) = LIF_StSol([2*k/Nfr;0],P,1.25);
    disp(k)
end
%%
k=23;plot([StFR0(:,k) StFR1(:,k) StFR4pi(:,k) StFR7(:,k)])
%%

plot(2*(1:Nfr)/Nfr,StFR0(end,:)*P.Fsteady(end)/P.Vres)
ylabel('Firing rate (spikes neruon^-^1 s^-^1)')
xlabel('Input excitatory conductance (a.u.)')



%% Single neuron Plots

plot(dt1:dt1:10*N*dt1,VL(:,1:6))
hold on
plot(dt:dt:N*dt, .05*LFP(2,:)-10,'k','LineWidth',2)
hold off
xlabel('\fontsize{24} Time (s)')
ylabel('\fontsize{24} Voltage (mV)')
ylim([-95 5])


