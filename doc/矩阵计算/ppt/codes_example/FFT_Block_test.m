% Block FFT
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

% testing kron(F,I)
x = rand(n) + i*rand(n);
y1 = kron(F,eye(n)) * x(:);
y2 = fft(x.').';
fprintf('kron(F,I): ||y1-y2||_inf=%.4e\n', norm(y1-y2(:),inf))

% testing kron(F'/n,I)
x = rand(n) + i*rand(n);
y1 = kron(F'/n,eye(n)) * x(:);
y2 = ifft(x.').';
fprintf('kron(F''/n,I): ||y1-y2||_inf=%.4e\n', norm(y1-y2(:),inf))

