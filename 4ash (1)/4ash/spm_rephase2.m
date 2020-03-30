function [y0] = spm_rephase2(y,x,n,f,Ang)
% k=2; I = M.opt.tu(:,k) > M.opt.tu(end,k) - 1/M.fu(k); y=squeeze(dxdp(:,:,:,k));x=M.opt.tu(:,k); f= M.fu(k); n =100;

% I = M.opt.tu(:,k) > M.opt.tu(end,k) - 1/M.fu(k); 

s     = size(y);
xi    = (x(end-5)-(1/f):1/f/n:x(end-5));
y     = reshape(y,[prod(s(1:end-1)) s(end)]);


yi = interp1(x,y',xi); %yi = [yi(301:end,:);yi(1:300,:)];
yi = yi - (0:1/n:1)'*(yi(end,:)-yi(1,:));
% disp([norm(yi(end,:)-yi(1,:)) norm(yi(2,:)-yi(1,:))])


yf  = fft(yi(2:end,:));

if mod(n,2)
   
    phaseshift = exp(1i*[0:(n-1)/2,-(n-1)/2:-1]*(angle(sum(yf(2,1)))-Ang))'; % for n odd
    phaseshift(ceil(1*n/6+1):floor(5*n/6)) = 0;
%     phaseshift(1) = 0;
else
    
    phaseshift = exp(1i*[0:n/2-1,-n/2:-1]*(angle(sum(yf(2,1)))-Ang))'; % for n even
    phaseshift(ceil(1*n/6):floor(5*n/6)) = 0;
%     phaseshift(1) = 0;
end

yf0 = bsxfun(@times,phaseshift,yf);

y0  = real(ifft(yf0));

y0 = reshape(y0',[s(1: end-1) n]);

end
