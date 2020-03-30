  function Y = collect(c,field)
        Y = [];
        for k = 1:numel(c)
            d = shiftdim(getfield(c{k},field),-1);
            Y = [Y;d];
            
        end
    end