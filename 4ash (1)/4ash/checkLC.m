function bol = checkLC(Currents,t,Kend,t0,r0,tol)

bol = 1;

if t(Kend)-t0 <1/200
    bol = 0;
    return
end

I = t>t0;

D = sum(bsxfun(@minus,Currents(:,Kend),Currents(:,I)).^2);

Tol = tol * sum(Currents(:,Kend).^2);

if max(D) < Tol
    bol=0;
    return
end

[pks,locs] = findpeaks(-D);

if length(pks) > length(D)/8
    bol = 0;
    return
end
   
end
