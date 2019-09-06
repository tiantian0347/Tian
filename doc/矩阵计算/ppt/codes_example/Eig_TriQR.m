clear all;
n = 5;
A = randn(n);
T = hess(A+A');  % 对称三对角化

% T = diag([1.263,-0.82812,-3.1883],-1) ...
%     + diag([0.24929,0.96880,0.48539,-0.91563]) ...
%     + diag([1.263,-0.82812,-3.1883],1);

format short e;

% Perform QR iteration with Wilkinson's shift
for k = 1 : 4
    % Compute the shift
    sub_matrix = T(n-1:n,n-1:n);
    Eig_sub_matrix = eig(sub_matrix);
    if ( abs(T(n,n) - Eig_sub_matrix(1)) < abs(T(n,n) - Eig_sub_matrix(2)) )
        shift = Eig_sub_matrix(1);
    else
        shift = Eig_sub_matrix(2);
    end
    
    % Perform QR iteration
    [Q,R] = qr(T-shift*eye(n));
    T = R*Q + shift*eye(n);
    % Enforce symmetry explicitly
    T = tril(triu(T,-1),1); T = (T+T')/2;
    
    %   fprintf('k=%d,T(n,n-1)=%.4e,T(n,n)=%.4e\n', k,T(n,n-1),T(n,n));
    fprintf('k = %d \n',k)
    T
end

