% ��λ�Ƶ� QR �����㷨��ʾ
clear all; 
pause on;
format short e;
iter_max = 100;

Lam = [9, 5, 3, 1];            % ����ֵ
n = length(Lam);
rng(2015)
X = rand(n);

A = (X*diag(Lam))/X             % �� Lam Ϊ����ֵ�ľ���
pause

tol = max(abs(A(:)))/1e6;       % �����ǲ�����С�� tol �Ľ�ֱ����Ϊ 0

% �����ж�
L = A - triu(A); % A �������ǲ���
if (max(abs(L(:))) < tol) 
    return;
end

% ȷ��λ�� sigma
idx = n; 
while (idx > 1)
    if ( sum(abs(L(idx,:)))>0 )
        sigma = A(idx,idx); break;
    else
        idx = idx - 1;
    end
end

% ��ʼ����
for k = 1 : iter_max
    [Q,R] = qr(A-sigma*eye(n));
    A = R*Q + sigma*eye(n);
    
    L = A - triu(A);            % A �������ǲ��� 
    L(find(abs(L)<tol)) = 0;    % ������ֵС�� tol ����Ϊ 0 
    A = L + triu(A);
    
    fprintf('k = %d, sigma=%.4e\n',k,sigma); 
    A    % ��� A_k
    pause 
    
    % �����ж�
    if (max(abs(L(:))) < tol)   
        break;
    end
    
    % ȷ��λ��
    while (idx>1)
        if (sum(abs(L(idx,:)))>0)
            sigma = A(idx,idx); break;
        end
        idx = idx - 1;
    end
end

