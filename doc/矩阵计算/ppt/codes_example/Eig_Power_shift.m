% ��λ�Ƶ��ݵ���������ʾ
clear all; 
itermax = 20;
Lam = [9, 5, 3, 1];   % ��������ֵ
n = length(Lam);
X = rand(n);
A = (X*diag(Lam))/X;  % �� Lam Ϊ����ֵ�ľ���
x0 = ones(n,1)/sqrt(n);  % ������ʼ����

% power method
x = x0;
for k = 1 : itermax
    x = A*x;
    x = x / norm(x);   % ������������
    err(k) = abs(dot(A*x,x) - 9);  % ��������ֵ�����
end

% shifted power method
x = x0;
sigma = 3;
B = A - sigma*eye(n);
for k = 1 : itermax
    x = B*x;
    x = x / norm(x);   % ������������
    err_s(k) = abs(dot(B*x,x) + sigma - 9);  % ��������ֵ�����
end
 
% ��ͼ�Ƚ�
tt = 1:itermax;
semilogy(tt,err,'ob', tt,err_s,'+r');


