% DST
clear;

n = 2^8;

% testing dst
x = rand(n,1);
y1 = zeros(n,1);
for k = 1:n
    y1(k) = sin((k*pi/(n+1))*[1:n]) * x;
end
y2 = dst(x);
fprintf(' DST: ||y1-y2||_inf=%.4e\n', norm(y1-y2,inf))

% testing idst
y1 = zeros(n,1);
for k = 1:n
    y1(k) = sin((k*pi/(n+1))*[1:n]) * x;
    y1(k) = y1(k) *(2/(n+1)); 
end
y2 = idst(x);
fprintf('IDST: ||y1-y2||_inf=%.4e\n', norm(y1-y2,inf))


% 
% % Fourier matrix
% F = zeros(n);
% w = cos(2*pi/n) - i*sin(2*pi/n);
% for k = 1:n
%     for j = 1:n
%         F(k,j) = w^((k-1)*(j-1));
%     end
% end