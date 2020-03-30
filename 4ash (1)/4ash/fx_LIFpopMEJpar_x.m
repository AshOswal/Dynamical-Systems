function [ dxdp, x ] = fx_LIFpopMEJpar_x(P0,M,Vp)

% cflag = 0;
M.opt.Drawplots = 0;
Vp1      = [Vp zeros(size(Vp,1),1)];
% M.opt.tu = zeros(NT0,size(Vp1,2));
fu       = zeros(1,size(Vp1,2));
cu       = zeros(1,size(Vp1,2));

s          = [size(M.Currents) size(Vp1,2)];
s(2)       = M.opt.Ncoefs;
LCrephased = zeros(s);
FR         = zeros(size(M.opt.popFiringRates,1), size(Vp1,2));

fprintf('\n Computing Gradients:')
ndp = size(Vp1,2);


parfor j =1:ndp
    
      Pj  = spm_ExpP(spm_unvec(spm_vec(P0) +  M.opt.dpsize*Vp1(:,j),P0));
      Mj  = spm_lifpopsys_LC_adapt(Pj,M);
      
      [~, ~, Mj] = fx_LIFpopMEJpar(Pj,Mj);
      
      IXi = Mj.J.IX(imag(Mj.J.S(Mj.J.IX))>.2*pi);
      lamb_max = Mj.J.S(IXi(end));
      
      fu(j)  = abs(imag(lamb_max))/(2*pi);
      cu(j)  = real(lamb_max)-1;% cu(j) = cu(j);%*(cu(j)<0);
      
      Mj.opt.N = 2;
      Mj = spm_DCM_lifpopsys_LC_int(Pj,Mj,0);

      LCrephased(:,:,j) = 0*(ones(Mj.opt.Ncoefs,1)*(Mj.Currents(:,2))')';
      FR(:,j) = Mj.GV(M.opt.popFiringRates,1,2);
      
end

x.lc =  0*permute(LCrephased(:,:,end),[2 1]); % low frequencies are removed from currents
x.f  =  fu(end)*ones(M.opt.Ncoefs,1);
x.fr =  ones(M.opt.Ncoefs,1)*FR(:,end)';
x.c  =  cu(end)*ones(M.opt.Ncoefs,1);

dxdp.lc = 0*permute(bsxfun(@minus,LCrephased(:,:,1:end-1),LCrephased(:,:,end)),[2 1 3]); % low frequencies are removed from currents
dxdp.f  = ones(M.opt.Ncoefs,1)*fu(1:end-1) -  fu(end);
dxdp.f  = reshape(dxdp.f,[size(dxdp.f,1) 1 size(dxdp.f,2)]);
dxdp.fr = bsxfun(@minus,FR(:,1:end-1),FR(:,end));
dxdp.fr = repmat(reshape(dxdp.fr,[1 size(dxdp.fr)]),[M.opt.Ncoefs 1 1]);
dxdp.c  = ones(M.opt.Ncoefs,1)*cu(1:end-1) -  cu(end);
dxdp.c  = reshape(dxdp.c,[size(dxdp.c,1) 1 size(dxdp.c,2)]);

dxdp    = cat(2,dxdp.lc,dxdp.f,dxdp.fr,dxdp.c)./M.opt.dpsize;



