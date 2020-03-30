function M = mate(M1,M2)

% M = M1;
% 
% C = M1.Cp + M2.Cp;
% 
% Ind = find(diag(C)~=0);
% 
% P = spm_vec(M1.Ep);
% 
% P1 = spm_vec(M1.Ep);
% P2 = spm_vec(M2.Ep);
% 
% P(Ind) = (C(Ind,Ind)/M1.Cp(Ind,Ind))*P1(Ind) + (C(Ind,Ind)/M2.Cp(Ind,Ind))*P2(Ind);
% 
% M.P = P;


M  = M1;
C1 = M1.F*spm_inv(M1.Cp);
C2 = M2.F*spm_inv(M2.Cp);
C  = (C1 + C2);

P = spm_vec(M1.Ep);

P1 = spm_vec(M1.Ep);
P2 = spm_vec(M2.Ep);

P = C\(C1*P1 + C2*P2);

M.P  = P;
M.Cp = C;


end