% Rayleigh 商迭代方法演示
clear all; 

itermax = 5;
Lam = [9, 5, 3, 1];   % 给定特征值
n = length(Lam);
rng(2017); 
X = rand(n);

A = (X*diag(Lam))/X;  % 以 Lam 为特征值的矩阵

x0 = ones(n,1)/sqrt(n);  % 迭代初始向量

% Rayleigh 商迭代
x = x0;
sigma = dot(A*x,x);
for k = 1 : itermax
    x = (A-sigma*eye(n)) \ x;
    x = x / norm(x);   
    sigma  = dot(A*x,x);
    [tmp,idx] = min(abs(Lam - sigma)); % 找出最近的特征值
    err(k) = abs( sigma - Lam(idx));   % 与最近特征值之间的误差
end
 
% 绘图
tt = 1:itermax;
semilogy(tt,err,'ob');


