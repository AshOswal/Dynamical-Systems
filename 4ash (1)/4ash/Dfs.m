function [f] = Dfs(f)

% Discrete Fourier series for real valued vectors

n = size(f,1);
f = fft(f)/n;


if mod(n,2)
    f(1:(n+1)/2,:) = real(f(1:(n+1)/2,:)); % for n odd
    f((n+1)/2+1:end,:) = imag(f((n+1)/2+1:end,:));
else
    f(1:(n)/2,:) = real(f(1:(n)/2,:)); % for n even
    f((n)/2+1:end,:) = imag(f((n)/2+1:end,:));
end

end