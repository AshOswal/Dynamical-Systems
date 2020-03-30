function fx_venous_outflow


P        = [];
P.T      = [0  6000];
P.Cpan   = 0.205;     % ml/mmHg
P.DCpa1  = 2.87;      % ml/mmHg
P.DCpa2  = 0.164;     % ml/mmHg
P.Gaut   = 3;
P.Ke     = 0.077;     % /ml
P.Kr     = 13.1*1e3;  % mlHg3.s/m
P.Kven   = 0.155;     % /ml
P.a      = 100;       % mmHg
P.icn    = 9.5;       % mmHg
P.Pa     = 58.9;      % mmHg
P.Pv     = 14.1;      % mmHg
P.Pv1    = -2.5;      % mmHg
P.Pvs     = 6;        % mmHg
P.Qn     = 12.5;      % ml/s
P.Ro     = 526.3*8;     % mmHg.s.ml-1
P.Rf     = 2.38*1e3;  % mmHg.s.ml-1
P.Rla    = 0.6;       % mmHg.s.ml-1
P.Rpv    = 0.880;     % mmHg.s.ml-1
P.Rvs1   = 0.366;     % mmHg.s.ml-1
P.taut   = 20;        % s
P.xaut   = 2.16*1e-4; 

% Capacities all in ml/mmHg
P.Cvs    = 0.5;       % venous sinus
P.Cjr3   = 1;         % right superior jugular
P.Cjl3   = 1;         % left superior jugular
P.Cjr2   = 2.5;       % right middle jugular
P.Cjl2   = 2.5;       % left middle jugular
P.Cc3    = 0.7;       % superior central collateral;
P.Cc2    = 1.4;       % middle central collateral
P.Csvc   = 20;        % SVC lower
P.Cazy   = 0.5;       % azygous vein
P.Cvv    = 0.5;       % vertebral 

% Other basal parameters
P.A      = 0.8;
P.kjr3   = 11.0;
P.kjr2   = 13.0;
P.kjr1   = 6.9;
P.kjl3   = 11.0;
P.kjl2   = 13.0;
P.kjl1   = 6.9;

% Venous pressures and flows
P.Pc3    = 6;         % mmHg --- collateral superior tract
P.Pjr3   = 5.85;      % mmHg --- right jugular (superior tract)
P.Pjl3   = 5.85;      % mmHg --- left jugular  (superior tract)
P.Pc2    = 5.85;      % mmHg --- collateral middle tract
P.Pvv    = 5.8;       % mmHg --- vertebral veins
P.Pjr2   = 5.7;       % mmHg --- right jugular (middle tract)
P.Pjl2   = 5.7;       % mmHg --- left jugular (middle tract)
P.Plv    = 5.6;       % mmHg --- lumbar vein
P.Pazy   = 5.5;       % mmHg --- azygous vein
P.Psvc1  = 5.4;       % mmHg --- SVC (superior tract)
P.Psvc   = 5.2;       % mmHg --- SVC (inferior tract)
P.Pcv    = 5.0;       % mmHg --- CVP (right atrium)

P.Qn     = 12.5;      % ml/s --- total cerebral blood flow
P.Qjr3   = 5.85;      % ml/s --- right jugular venous flow
P.Qjl3   = 5.85;      % ml/s --- left jugular venous flow
P.Qvvr   = 0.40;      % ml/s --- right vertebral venous flow
P.Qvvl   = 0.40;      % ml/s --- left vertebral venous flow
P.Qex    = 5.00;      % ml/s --- external duct
P.Qcjr3  = 1.00;      % ml/s --- anastomotic connnection
P.Qcjr2  = 1.00;      % ml/s --- anastomotic connection
P.Qcjl3  = 1.00;      % ml/s --- anastomotic connection
P.Qcjl2  = 1.00;      % ml/s --- anastomotic connection
P.Qc2    = 3.00;      % ml/s --- middle collateral
P.Qazy1  = 0.40;      % ml/s --- azygous vein
P.Qlv    = 0.13;      % ml/s --- lumbar vein
P.Qrv    = 0.27;      % ml/s --- renal vein

% Conductances (G) inverse of resistance = flow/pressure gradient
P.Gjr3   =  P.Qjr3/(P.Pvs - P.Pjr3);                                        % mmHg/ml/s --- right jugular superior tract
P.Gjl3   =  P.Qjl3/(P.Pvs - P.Pjl3);                                        % mmHg/ml/s --- left jugular superior tract
P.Gc3    =  0;%0/(P.Pvs - P.Pc3);                                              % mmHg/ml/s --- superior collateral 
P.Gcjr3  =  P.Qcjr3/(P.Pc3 - P.Pjr3);                                       % mmHg/ml/s --- upper anastomotic segment
P.Gcjl3  =  P.Qcjl3/(P.Pc3 - P.Pjl3);                                       % mmHg/ml/s --- upper anastomotic segment
P.Gex    =  P.Qex/(P.Pa -  P.Pc3);                                          % mmHg/ml/s --- of external tract 
P.Gvvl   =  P.Qvvl/(P.Pvs - P.Pvv);                                         % mmHg/ml/s --- left vertebral
P.Gvvr   =  P.Qvvr/(P.Pvs - P.Pvv);                                         % mmHg/ml/s --- right vertebral
P.Gc2    =  P.Qc2/(P.Pc3 - P.Pc2);                                          % mmHg/ml/s --- middle collateral
P.Gcjr2  =  P.Qcjr2/(P.Pc2 - P.Pjr2);                                       % mmHg/ml/s --- middle anastomotic segment
P.Gcjl2  =  P.Qcjl2/(P.Pc2 - P.Pjl2); % mmHg/ml/s --- middle anastomotic segment

P.Qjr2   = (P.Qjr3 + P.Qcjr3);                                               % summation of flows to create Qjr2
P.Qjl2   = (P.Qjl3 + P.Qcjl3); % summation of flows to create Qjl2

P.Gjr2   = P.Qjr2/(P.Pjr3 - P.Pjr2);
P.Gjl2   = P.Qjl2/(P.Pjl3 - P.Pjl2);
P.Gc1    = (P.Qc2 - P.Qcjr2 - P.Qcjl2)/(P.Pc2 -  P.Pcv);
P.Gjr1   = (P.Qjr2 + P.Qcjr2)/(P.Pjr2 - P.Psvc1);
P.Gjl1   = (P.Qjl2 + P.Qcjl2)/(P.Pjl2 - P.Psvc1);
P.Gsvc1  = (P.Qjr2 + P.Qcjr2 + P.Qjl2 + P.Qcjl2)/(P.Psvc1  - P.Psvc);
P.Gsvc2  = (P.Qjr2 + P.Qcjr2 + P.Qjl2 + P.Qcjl2 + P.Qazy1)/(P.Psvc  - P.Pcv);
P.Gazy1  = P.Qazy1/(P.Pvv   - P.Pazy);
P.Gazy2  = (P.Qazy1 + P.Qlv)/(P.Pazy  - P.Psvc);
P.Gvv2   = (P.Qlv + P.Qrv)/(P.Pvv   - P.Plv);
P.Glv    = P.Qlv/(P.Plv   - P.Pazy); 


F        = [];
F.xaut   = P.xaut;
F.icn    = P.icn;
F.Pv     = P.Pv;
F.Pa     = P.Pa;
F.Pvs    = P.Pvs;
F.Pjr3   = P.Pjr3;
F.Pjr2   = P.Pjr2;
F.Pjl3   = P.Pjl3;
F.Pjl2   = P.Pjl2;
F.Pc3    = P.Pc3;
F.Pc2    = P.Pc2;
F.Psvc   = P.Psvc;
F.Pvv    = P.Pvv;
F.Pazy   = P.Pazy;

%options = odeset('RelTol', 1e-5,'AbsTol',1e-5);
[T, f]  = ode45(@(t,F) fx_icv(t,F,P),P.T,[P.xaut, P.icn, P.Pv, P.Pa, P.Pvs, P.Pjr3, P.Pjr2, P.Pjl3, P.Pjl2, P.Pc3, P.Pc2, P.Psvc, P.Pvv, P.Pazy ]);
figure;plot(T,f(:,1:5));legend({'xaut','ICP','PV','PA','PVs','PJR3','PJR2','PJL3','PJL2','PC3','PC2','PSVC','PVV','PAZY'});
% equations of motion dfdt

function [dfdt] = fx_icv(t,F,P)

% F - values of time dependent variables

% dfdt     = [dXautdt, dPicdt, dPvdt, dPadt, dPvsdt, dPjr3dt, dPjr2dt, dPjl3dt, dPjl2dt, dPc3dt, dPc2dt, dPsvcdt, dPvvdt, dPazydt];

% F.xaut - autoregulation state variable
% F.Pic  - intracranial pressure        
% F.Pv   - terminal ic venous pressure
% F.Pa   - pial arterioles
% F.Pvs  - venous sinus pressure
% F.Pjr3 - right jugular
% F.Pjr2 - right jugular (M)
% F.Pjl3 - left jugular
% F.Pjl2 - left jugular  (M)
% F.Pc3  - superior collateral
% F.Pc2  - middle collateral
% F.Psvc - SVC
% F.Pvv  - vertebral vein

xaut     = F(1);
Pic      = F(2);
Pv       = F(3);
Pa       = F(4);
Pvs      = F(5);
Pjr3     = F(6);
Pjr2     = F(7);
Pjl3     = F(8);
Pjl2     = F(9);
Pc3      = F(10);
Pc2      = F(11);
Psvc     = F(12);
Pvv      = F(13);
Pazy     = F(14);


% autoregulation phase
if xaut < 0
   DCpa = P.DCpa1;
else
   DCpa = P.DCpa2; 
end
kCpa = DCpa/4;

% compliance & resistance of pial arterioles

Cpa = ((P.Cpan - DCpa/2) + (P.Cpan + DCpa/2)*exp(-xaut/kCpa)) / (1 + exp(-xaut/kCpa));
Rpa = (P.Kr*(P.Cpan^2))/(((Pa - Pic)*Cpa)^2);


% Pcap - capillary pressure
Pcap = ((Pv/P.Rpv)+(Pa/(Rpa/2))+(Pic/P.Rf)) / ((1/P.Rpv)+(1/(Rpa/2))+(1/P.Rf));

% Q    - cerebral blood flow
Q    = (Pa - Pcap)/(Rpa/2);

% flow dxdt
dXautdt = (1/P.taut)*(-xaut + P.Gaut*((Q - P.Qn)/P.Qn));

% dCpadt - change in arterial compliance -- analytically derived from Cpa
dCpadt  = ((DCpa/kCpa) * exp(-xaut/kCpa) * -dXautdt) / ((1 + exp(-xaut/kCpa))^2);

% CSF formation and absorption
if Pcap > Pic
   Qf = (Pcap - Pic)/P.Rf;
else
   Qf = 0;
end

if Pic > Pvs
   Qo = (Pic - Pvs)/P.Ro;
else
   Qo = 0;
end

% venous sinus resistance
if Pv - Pvs > 0 
   Rvs = ((Pv - Pvs)/(Pv - Pic))*P.Rvs1;
else
   Rvs = P.Rvs1;
end

% Cvi  - compliance of intracranial veins
Cvi  = 1/(P.Kven*(Pv - P.Pv1 - Pic));

% Cic  - intracranial compliance
Cic  = 1/(P.Ke*Pic); 

% dPicdt - computation

% d(Pv-Pic)/dt
dPvIcdt = (1/Cvi)*(((Pcap - Pv)/P.Rpv)  - ((Pv - Pvs)/Rvs));


% d(Pa-Pic)/dt
dPaIcdt = (1/Cpa)* ((P.a - Pa)/(P.Rla + Rpa/2) - (Pa - Pcap)/(Rpa/2) - dCpadt*(Pa - Pic));

dPicdt  = (1/Cic)* (dPaIcdt*Cpa + dPvIcdt*Cvi +    dCpadt*(Pa - Pic) + Qf - Qo );

% dPv/dt
dPvdt   = dPvIcdt + dPicdt;

% dPa/dt
dPadt   = dPaIcdt + dPicdt;

%% TREAT THIS SEPARATE - NO UPDATING
%  computation of conductances from physiological data

P.Gvs    = 1/Rvs;   

    % mmHg/ml/s --- termnal intracranial veins
P.Go     = 1/P.Ro;                                                         % mmHg/ml/s --- to CSF outflow

% update conductances of jugular system in supine condition
P.Gjr3   = P.kjr3* ((1 + (2/pi)*atan((Pvs - 0)/P.A)))^2;
P.Gjl3   = P.kjl3* ((1 + (2/pi)*atan((Pvs - 0)/P.A)))^2;
P.Gjr2   = P.kjr2* ((1 + (2/pi)*atan((Pjr3 - 0)/P.A)))^2;
P.Gjl2   = P.kjl2* ((1 + (2/pi)*atan((Pjl3 - 0)/P.A)))^2;
P.Gjr1   = P.kjr1* ((1 + (2/pi)*atan((Pjr2 - -6.5)/P.A)))^2;
P.Gjl1   = P.kjl1* ((1 + (2/pi)*atan((Pjl2 - -6.5)/P.A)))^2;

% if t >=2
%    P.Gjr3 = 0.01*P.Gjr3;
%    P.Gjl3 = 0.01*P.Gjl3;
% end
%% Venous outflow equations of motion dfdt

dPvsdt   = 1/P.Cvs * ((Pv - Pvs)*P.Gvs - (Pvs - Pic)*P.Go - (Pvs - Pjr3)*P.Gjr3 - (Pvs - Pjl3)*P.Gjl3 - (Pvs - Pc3)*P.Gc3 - (Pvs - Pvv)*P.Gvvl - (Pvs - Pvv)*P.Gvvr);

dPjr3dt  = 1/P.Cjr3* ((Pvs - Pjr3)*P.Gjr3 - (Pjr3 - Pc3)*P.Gcjr3 - (Pjr3 - Pjr2)*P.Gjr2);

dPjr2dt  = 1/P.Cjr2* ((Pjr3 - Pjr2)*P.Gjr2 - (Pjr2 - Pc2)*P.Gcjr2 - (Pjr2 - P.Psvc1)*P.Gjr1);

dPjl3dt  = 1/P.Cjl3* ((Pvs - Pjl3)*P.Gjl3 - (Pjl3 - Pc3)*P.Gcjl3 - (Pjl3 - Pjl2)*P.Gjl2);

dPjl2dt  = 1/P.Cjl2* ((Pjl3 - Pjl2)*P.Gjl2 - (Pjl2 - Pc2)*P.Gcjl2 - (Pjl2 - P.Psvc1)*P.Gjl1);

dPc3dt   = 1/P.Cc3 * ((Pvs - Pc3)*P.Gc3  + (Pjr3 - Pc3)*P.Gcjr3 + (Pjl3 - Pc3)*P.Gcjl3   - (Pc3 - Pc2)*P.Gc2 +(P.a - Pc3)*P.Gex); % 

dPc2dt   = 1/P.Cc2 * ((Pc3 - Pc2)*P.Gc2  + (Pjr2 - Pc2)*P.Gcjr2 + (Pjl2 - Pc2)*P.Gcjl2 - (Pc2 - P.Pcv)*P.Gc1);

dPsvcdt  = 1/P.Csvc* ((P.Psvc1 - Psvc)*P.Gsvc1 + (Pazy - Psvc)*P.Gazy2 -(Psvc - P.Pcv)*P.Gsvc2);

dPvvdt   = 1/P.Cvv*  ((Pvs - Pvv)*P.Gvvl + (Pvs - Pvv)*P.Gvvr - (Pvv - Pazy)*P.Gazy1 - (Pvv - P.Plv)*P.Gvv2);

dPazydt  = 1/P.Cazy* ((Pvv - Pazy)*P.Gazy1 + (P.Plv - Pazy)*P.Glv - (Pazy - Psvc)*P.Gazy2);


dfdt     = [dXautdt, dPicdt, dPvdt, dPadt, dPvsdt, dPjr3dt, dPjr2dt, dPjl3dt, dPjl2dt, dPc3dt, dPc2dt, dPsvcdt, dPvvdt, dPazydt]';



