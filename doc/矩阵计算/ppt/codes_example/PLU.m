function [A,p] = PLU(A)
%
% 列主元 LU 分解
% PLU : A(p,:) = L * U
%

[n,n] = size(A);
p = 1 : n; 
for k = 1 : n-1
    [a_max,l] = max(abs(A(k:n,k)));
    if a_max == 0
        fprintf('Error: 第 %d 步的列主元为 0!\n', k); return;
    end
    l = l+k-1;
    if l~=k
        t = A(k,:); A(k,:) = A(l,:); A(l,:) = t;
        tmp = p(k); p(k) = p(l); p(l) = tmp;
    end
    A(k+1:n,k) = A(k+1:n,k)/A(k,k);
    A(k+1:n,k+1:n) = A(k+1:n,k+1:n) - A(k+1:n,k)*A(k,k+1:n);
end