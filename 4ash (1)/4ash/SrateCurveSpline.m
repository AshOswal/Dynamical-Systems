function pp = SrateCurveSpline(P,In,res,Emax)


GoV = -1:res:Emax;
SrateV = zeros(length(GoV),1);
VrestV = SrateV;
for k=1:length(GoV)
    
    Go = [GoV(k) In];
    [Qdpdt] = fx_LIFpopME(Go,P);
    Qdpdt = Qdpdt(1:end-length(Go),1:end-length(Go));
    [V S] = eig(full(Qdpdt));
    [B,IX] = sort(real(diag(S)));
%     subplot(2,1,1)
%     plot(B)
%     subplot(2,1,2)
%     plot(P.VV,abs([V(1:P.LVV,IX(end))]))
%     input('')
    Vss       = V(:,IX(end));
    Srate     = Vss(end)/sum(Vss)*P.Fsteady(end);
    SrateV(k) = Srate;
end

plot(GoV,SrateV)
ylabel('Firing Rate (#spikes/neuron/s)','FontSize',14)
% pp = interp1(VrestV,SrateV,'spline','pp');

end