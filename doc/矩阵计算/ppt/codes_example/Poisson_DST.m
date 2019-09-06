%
% DST for the 2-D Poisson problem
%

clear all;
close all;
N = 2^4;
h = 1/(N+1); % step size

% initial guess + boundary condition
v0 = zeros(N+2,N+2); 
x = 0:h:1;  y = 0:h:1;
v0(1,:) = 0.25 * (y .* y);
v0(N+2,:) = 0.25 * (1 + y .* y);
v0(:,1) = 0.25 * (x .* x);
v0(:,N+2) = 0.25 * (x .* x + 1);

[X,Y] = meshgrid(x);
vt = 0.25*(X.*X + Y.*Y); % true solution
VT = vt(2:end-1,2:end-1); % interior points

% eigenvalues of 1-D and 2-D Poisson equation
Lam1D = 2-2*cos((1:N)*(pi/(N+1)));
Lam1D = Lam1D(:);
Lam2D = kron(ones(N,1),Lam1D) + kron(Lam1D,ones(N,1));
Lam2D = reshape(Lam2D,N,N);

% right hand side
F = (-1*h*h)*ones(N); % h^2*f(x,y) 
F(1,:) = F(1,:) + vt(1,2:end-1);
F(:,1) = F(:,1) + vt(2:end-1,1);
F(:,end) = F(:,end) + vt(2:end-1,end);
F(end,:) = F(end,:) + vt(end,2:end-1);

% solve 2-D poisson equation with DST
% U = idst(F).';
% U = idst(U).';
% U = U ./ reshape(Lam2D,N,N);
% U = dst(U).';
% U = dst(U).';
U = dst(dst((idst(idst(F).').') ./ Lam2D).').';
fprintf('DST-2D: error=%.4e\n', norm(U(:)-VT(:),inf));

% poicalc -- fast solver for Poisson equation (Matlab function)
U2 = poicalc(F(:),h,h,N,N);
fprintf('poicalc: error=%.4e\n', norm(U(:)-U2(:),inf));

% direct method
if (N<=2^6)  % the problem can not be too large
    T1D = 2*eye(N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1);
    T2D = kron(T1D,eye(N)) + kron(eye(N),T1D);
    UT = T2D\F(:);
    fprintf('direct: error=%.4e\n', norm(U(:)-UT,inf));
end
