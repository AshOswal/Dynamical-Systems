function fx_iih_venous


P        = [];
P.T      = [0 :1 :3500];
P.Cpan   = 0.205;     % ml/mmHg
P.DCpa1  = 2.87;      % ml/mmHg
P.DCpa2  = 0.164;     % ml/mmHg
P.Gaut   = 3;
P.Ke     = 0.077;     % /ml
P.Kr     = 13.1*1e3;  % mlHg3.s/ml
P.Kven   = 0.155;     % /ml
P.a      = 100;       % mmHg
P.icn    = 9.5;       % mmHg
P.Pa     = 58.9;      % mmHg
P.Pv     = 14.1;      % mmHg
P.Pv1    = -2.5;      % mmHg
P.Pvs     = 6.0;        % mmHg
P.Qn     = 12.5;      % ml/s
P.Ro     = 4*526.3;     % mmHg.s.ml-1
P.Rf     = 2.38*1e3;  % mmHg.s.ml-1
P.Rla    = 0.6;       % mmHg.s.ml-1
P.Rpv    = 0.880;     % mmHg.s.ml-1
P.Rvs1   = 0.366;     % mmHg.s.ml-1
P.taut   = 20;        % s
P.xaut   = 2.16*1e-4; 

% Capacities all in ml/mmHg
P.Cvs    = 0.5;       % venous sinus
P.Ctr    = 0.5;       % right transverse sinus
P.Ctl    = 0.5;       % left  transverse sinus
P.Cjr    = 1;       % right jugular
P.Cjl    = 1;       % left jugular


% Venous pressures and flows
P.Ptr    = 5.9;       % mmHg --- right transverse sinus
P.Ptl    = 5.9;       % mmHg --- left transverse sinus 
P.Pjr    = 5.85;       % mmHg --- right jugular 
P.Pjl    = 5.85;       % mmHg --- left jugular 

P.Qn     = 12.5;        % ml/s --- total cerebral blood flow
P.Qtr    = P.Qn/2;      % ml/s --- right venous flow
P.Qtl    = P.Qn - P.Qtr;% ml/s --- left ar venous flow

P.Qjr    =  P.Qtr;                                             
P.Qjl    =  P.Qtl; 

F        = [];
F.xaut   = P.xaut;
F.icn    = P.icn;
F.Pv     = P.Pv;
F.Pa     = P.Pa;
F.Pvs    = P.Pvs;
F.Ptr   =  P.Ptr;
F.Ptl   =  P.Ptl;

P.A      = 0.8;
P.kjr3   = 11.0;
P.kjr2   = 13.0;
P.kjl3   = 11.0;
P.kjl2   = 13.0;

P.Gtr   = P.Qtr/4*(P.Pvs - P.Ptr);                                   % mmHg/ml/s --- right jugular superior tract
P.Gtl   = P.Qtl/1.2*(P.Pvs - P.Ptl); 


P.Gjr   = P.kjr3* ((1 + (2/pi)*atan((P.Ptr - 0)/P.A)))^2;
P.Gjl   = P.kjl3* ((1 + (2/pi)*atan((P.Ptl- 0)/P.A)))^2;

options = odeset('RelTol', 1e-4,'AbsTol',1e-4);
[T, f]  = ode45(@(t,F) fx_icv(t,F,P),P.T,[P.xaut, P.icn, P.Pv, P.Pa, P.Pvs, P.Ptr, P.Ptl],options);
figure;plot(T,f);legend({'xaut','ICP','Venous pressure','Pial arterial pressure','Venous sinus pressure','Right transverse sinus pressure','Left transverse sinus pressure'});
% equations of motion dfdt


function [dfdt] = fx_icv(t,F,P)

% F - values of time dependent variables

% dfdt     = [dXautdt, dPicdt, dPvdt, dPadt, dPvsdt, dPtr, dPtl, dPjr, dPjl];

% F.xaut - autoregulation state variable
% F.Pic  - intracranial pressure        
% F.Pv   - terminal ic venous pressure
% F.Pa   - pial arterioles
% F.Pvs  - venous sinus pressure
% F.Pjr3 - right jugular
% F.Pjl3 - left jugular

xaut     = F(1);
Pic      = F(2);
Pv       = F(3);
Pa       = F(4);
Pvs      = F(5);
Ptr      = F(6);
Ptl      = F(7);

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
P.Go     = 1/P.Ro;       % mmHg/ml/s --- to CSF outflow

                                      % mmHg/ml/s --- left jugular superior tract
 

%% Venous outflow equations of motion dfdt

dPvsdt   = 1/P.Cvs * ((Pv - Pvs)*P.Gvs - (Pvs - Pic)*P.Go - (Pvs - Ptr)*P.Gtr - (Pvs - Ptl)*P.Gtl);

dPtrdt  = 1/P.Cjr* ((Pvs - Ptr)*P.Gtr  - (Ptr - P.Pjr)*P.Gjr);

dPtldt  = 1/P.Cjl* ((Pvs - Ptl)*P.Gtl  - (Ptl - P.Pjr)*P.Gjl);

dfdt     = [dXautdt, dPicdt, dPvdt, dPadt, dPvsdt, dPtrdt, dPtldt]';
if any(isnan(dfdt)) 
   keyboard;
end
if any(isinf(dfdt)) 
clear    keyboard;
end
