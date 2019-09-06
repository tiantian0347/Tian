% QR 迭代算法演示
clear all; 
pause on;
iter_max = 100;

A=[ 11.2  12.6    1.0  -13.7
    4.3    6.6    1.5   -6.4
    8.3   17.5    4.5  -17.6
    6.1    3.0    1.7   -4.3 ]

tol = max(abs(A(:)))/1e6;       % 下三角部分中小于 tol 的将直接设为 0

for k = 1 : iter_max
    [Q,R] = qr(A);
    A = R*Q;
    
    L = A - triu(A);            % 取 A 的下三角部分 
    L(find(abs(L)<tol)) = 0;    % 将绝对值小于 tol 的设为 0 
    A = L + triu(A);
    
    fprintf('k = %d\n',k); A    % 输出 A_k
    pause 
    
    if (max(abs(L(:))) < tol)   % 收敛判断
        break;
    end
end

