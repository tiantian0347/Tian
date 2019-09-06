% 二维离散 Poisson 方程的特征值
%
clear all;
close all;
tol = 1e-16;
N = 6 ;
T1D = 2*eye(N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1);
A = kron(eye(N),T1D) + kron(T1D,eye(N));

% 计算出所有特征值
Eig_A = eig(A);

% Lanczos 方法
n = N*N;
m = n/2;  % m Krylov 子空间最大允许维数

figure(11)
hold on
axis([0,m+2,0,8]);
plot((m+1)*ones(n,1),Eig_A,'+');
tit = sprintf('n=%d, m=%d',n,m);
title(tit)

% 分配空间
V = zeros(n,m);
alpha = zeros(m,1);
beta  = zeros(m,1);

% 开始迭代
V(:,1) = rand(n,1); V(:,1) = V(:,1)/norm(V(:,1));
for j = 1 : m
    if j > 1
        w = A * V(:,j) - beta(j-1)*V(:,j-1);
    else
        w = A * V(:,j);
    end
    alpha(j) = w.' * V(:,j);
    w = w - alpha(j) * V(:,j);
    beta(j) = norm(w);
    if beta(j) < tol, break, end
    V(:,j+1) = w / beta(j);
    T = diag(alpha(1:j)) + diag(beta(1:j-1),1) + diag(beta(1:j-1),-1);
    Eig_T = eig(T);
    plot(j*ones(j,1),Eig_T,'+');
    plot([j,j],[min(Eig_T),max(Eig_T)],'or');
    pause
end


