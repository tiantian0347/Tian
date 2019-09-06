% QR �����㷨��ʾ
clear all; 
pause on;
format short e;
iter_max = 100;

Lam = [9, 5, 3, 1];            % ����ֵ
n = length(Lam);
rng(2015); 
X = rand(n);

A = (X*diag(Lam))/X             % �� Lam Ϊ����ֵ�ľ���
pause

tol = max(abs(A(:)))/1e6;       % �����ǲ�����С�� tol �Ľ�ֱ����Ϊ 0

for k = 1 : iter_max
    [Q,R] = qr(A);
    A = R*Q;
    
    L = A - triu(A);            % A �������ǲ��� 
    L(find(abs(L)<tol)) = 0;    % ������ֵС�� tol ����Ϊ 0 
    A = L + triu(A);
    
    fprintf('k = %d\n',k); A    % ��� A_k
    pause 
    
    % �����ж�
    if (max(abs(L(:))) < tol)   
        break;
    end
end

