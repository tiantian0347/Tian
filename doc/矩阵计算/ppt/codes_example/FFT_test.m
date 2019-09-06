% FFT
clear;

n = 2^6;

% Fourier matrix
F = zeros(n);
w = cos(2*pi/n) - i*sin(2*pi/n);
for k = 1:n
    for j = 1:n
        F(k,j) = w^((k-1)*(j-1));
    end
end

% testing fft
x = rand(n,1)+i*rand(n,1);
y1 = F*x;
y2 = fft(x);
fprintf(' FFT: ||y1-y2||_inf=%.4e\n', norm(y1-y2,inf))

% testing ifft
y1 = (F' * x)/n;
y2 = ifft(x);
fprintf('IFFT: ||y1-y2||_inf=%.4e\n', norm(y1-y2,inf))


