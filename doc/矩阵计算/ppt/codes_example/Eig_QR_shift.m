% 带位移的 QR 迭代算法演示
clear all; 
pause on;
format short e;
iter_max = 100;

Lam = [9, 5, 3, 1];            % 特征值
n = length(Lam);
rng(2015)
X = rand(n);

A = (X*diag(Lam))/X             % 以 Lam 为特征值的矩阵
pause

tol = max(abs(A(:)))/1e6;       % 下三角部分中小于 tol 的将直接设为 0

% 收敛判断
L = A - triu(A); % A 的下三角部分
if (max(abs(L(:))) < tol) 
    return;
end

% 确定位移 sigma
idx = n; 
while (idx > 1)
    if ( sum(abs(L(idx,:)))>0 )
        sigma = A(idx,idx); break;
    else
        idx = idx - 1;
    end
end

% 开始迭代
for k = 1 : iter_max
    [Q,R] = qr(A-sigma*eye(n));
    A = R*Q + sigma*eye(n);
    
    L = A - triu(A);            % A 的下三角部分 
    L(find(abs(L)<tol)) = 0;    % 将绝对值小于 tol 的设为 0 
    A = L + triu(A);
    
    fprintf('k = %d, sigma=%.4e\n',k,sigma); 
    A    % 输出 A_k
    pause 
    
    % 收敛判断
    if (max(abs(L(:))) < tol)   
        break;
    end
    
    % 确定位移
    while (idx>1)
        if (sum(abs(L(idx,:)))>0)
            sigma = A(idx,idx); break;
        end
        idx = idx - 1;
    end
end

