function out=Nest(varargin)
func=varargin{1};
out=varargin{3};
for k=1:varargin{2}
    out=func(out);
end
end