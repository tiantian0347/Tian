% 一维 Poisson 方程离散矩阵的特征向量
%
clear all;
N = 64;
T = 2*eye(N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1);
Lam = diag(2-2*cos((1:N)*(pi/(N+1))));
V = zeros(N); % 特征向量
for k = 1 : N
    V(:,k) = sin((1:N)'*(k*pi/(N+1)));
end

norm(T*V-V*Lam,'fro')
norm(V(:,2))^2 *(2/(N+1))

x=rand(N,1);
y1=2/(N+1)*(V*x);
y2=idst(x);
norm(y1-y2,'fro')
