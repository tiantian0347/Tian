% QR �����㷨��ʾ
clear all; 
pause on;
iter_max = 100;

A=[ 11.2  12.6    1.0  -13.7
    4.3    6.6    1.5   -6.4
    8.3   17.5    4.5  -17.6
    6.1    3.0    1.7   -4.3 ]

tol = max(abs(A(:)))/1e6;       % �����ǲ�����С�� tol �Ľ�ֱ����Ϊ 0

for k = 1 : iter_max
    [Q,R] = qr(A);
    A = R*Q;
    
    L = A - triu(A);            % ȡ A �������ǲ��� 
    L(find(abs(L)<tol)) = 0;    % ������ֵС�� tol ����Ϊ 0 
    A = L + triu(A);
    
    fprintf('k = %d\n',k); A    % ��� A_k
    pause 
    
    if (max(abs(L(:))) < tol)   % �����ж�
        break;
    end
end

