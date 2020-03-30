function Mout=NormalizeM(Min)
Mout = Min - mean(Min,2)*ones(1,size(Min,2));
Mout = Mout./(sqrt(sum(Min.^2,2))*ones(1,size(Min,2)));
end