% 带位移的幂迭代方法演示
clear all; 
itermax = 20;
Lam = [9, 5, 3, 1];   % 给定特征值
n = length(Lam);
X = rand(n);
A = (X*diag(Lam))/X;  % 以 Lam 为特征值的矩阵
x0 = ones(n,1)/sqrt(n);  % 迭代初始向量

% power method
x = x0;
for k = 1 : itermax
    x = A*x;
    x = x / norm(x);   % 近似特征向量
    err(k) = abs(dot(A*x,x) - 9);  % 近似特征值的误差
end

% shifted power method
x = x0;
sigma = 3;
B = A - sigma*eye(n);
for k = 1 : itermax
    x = B*x;
    x = x / norm(x);   % 近似特征向量
    err_s(k) = abs(dot(B*x,x) + sigma - 9);  % 近似特征值的误差
end
 
% 绘图比较
tt = 1:itermax;
semilogy(tt,err,'ob', tt,err_s,'+r');


