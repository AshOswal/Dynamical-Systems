function [r, to]= convergence_checkCSVD(GV,t,k,T0)
% GV = M.GV; t = M.opt.t; T0=M.opt.T0;
r=Inf;
to = t(k)-T0;

GV(:,k+1:end) = [];
t(k+1:end) = [];

GV(:,t<t(end)-T0)=[];
t(t<t(end)-T0)=[];


s = size(GV);
VecGV = reshape(GV,[],s(end));

D = sum(bsxfun(@minus,VecGV(:,end),VecGV(:,1:end)).^2);

% [pks,locs] = findpeaks(-D);

[C, I] =  min(D(1:end-10));%max(pks);

% delta  = max(abs(D(locs(I)) - D(locs(I)-1:locs(I)+1)));
% delta2 = abs(D(end)-D(end-1));
% locs(pks<-delta+C) = [];
% pks(pks<-delta+C) = [];
% 
% locs(pks<-delta2+C) = [];


try
%     I = locs(end);
    
    curvexy = VecGV(:,I-3:I+3)';
    mapxy   = VecGV(:,end)';
    [xy,d,t_a]=distance2curve(curvexy,mapxy,'spline');

    
    r = norm(d);
    if nargout>1
        to = interp1(t(I-3:I+3),6*t_a+1,'spline');
    end
end




end



