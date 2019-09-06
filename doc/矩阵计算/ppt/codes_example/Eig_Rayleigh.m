% Rayleigh �̵���������ʾ
clear all; 

itermax = 5;
Lam = [9, 5, 3, 1];   % ��������ֵ
n = length(Lam);
rng(2017); 
X = rand(n);

A = (X*diag(Lam))/X;  % �� Lam Ϊ����ֵ�ľ���

x0 = ones(n,1)/sqrt(n);  % ������ʼ����

% Rayleigh �̵���
x = x0;
sigma = dot(A*x,x);
for k = 1 : itermax
    x = (A-sigma*eye(n)) \ x;
    x = x / norm(x);   
    sigma  = dot(A*x,x);
    [tmp,idx] = min(abs(Lam - sigma)); % �ҳ����������ֵ
    err(k) = abs( sigma - Lam(idx));   % ���������ֵ֮������
end
 
% ��ͼ
tt = 1:itermax;
semilogy(tt,err,'ob');


