function [f] = iDfs(f)

% inverse Discrete Fourier series for real valued vectors

n  = size(f,1);

if mod(n,2)
    
    f(2:(n+1)/2,:) = f(2:(n+1)/2,:) + 1i*f(end:-1:(n+1)/2+1,:); % for n odd
    f((n+1)/2+1:end,:) = conj(f((n+1)/2:-1:2,:));
      
else

    f(2:(n)/2,:) = f(2:(n)/2,:) + 1i*f(end:-1:(n)/2+2,:); % for n even
    f((n)/2+2:end,:) = conj(f((n)/2:-1:2,:));
        
end

f = real(fft(f))/n;

end