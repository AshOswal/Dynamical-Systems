function [y] = DawsonInt(x,f)

if f == 0
    y  = mfun('dawson',x);
else
    k = (1-f/2);
%     k = -((-2 + f)/(2*sqrt(f)));
    
    y1 = (1 - exp(-f*x.^2))./(2*x);
    y1(isnan(y1)) = 0; % prevent 0/0 at the origin
    
    y2 = exp(-x.^2).*sqrt(exp(x.^2) - 1).*sign(x)*k;
    y2(isnan(y2)|isinf(y2)) = 0; % prevent 0*Inf for large x
    
%     y1 = (1 - exp(-f*x.^2))./(2*x);
%     y1(isnan(y1)) = 0; % prevent 0/0 at the origin
%     
%     y2 = exp(-x.^2).*sqrt(exp(f*x.^2) - 1).*sign(x)*k;
%     y2(isnan(y2)|isinf(y2)) = 0; % prevent 0*Inf for large x
    
    y = y1 + y2;
end

end