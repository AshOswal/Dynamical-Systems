function [dVdt dGdt] = fx_LIF(V,Go,Gi,U,P)
dVdt = (sum((Go + P.DGL'.*P.DL').*(P.Vg'-V)) + P.gl*(P.Vl-V))/P.C;
if nargout > 1
    dGdt = (U-Gi)./P.T';
end
end