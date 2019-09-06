% 测试 SSOR 对参数 omega 的敏感性
%
%

clear all;
close all;
N = 8;
h = 1/(N+1); % step size

% initial guess
v0 = zeros(N+2,N+2);

% boundary condition
x = 0:h:1;  y = 0:h:1;
v0(1,:) = 0.25 * (y .* y);
v0(N+2,:) = 0.25 * (1 + y .* y);
v0(:,1) = 0.25 * (x .* x);
v0(:,N+2) = 0.25 * (x .* x + 1);

f = -1; % right-hand side f(x,y)

% true solution
[X,Y] = meshgrid(x);
vt = 0.25*(X.*X + Y.*Y);
norm_vt = norm(vt(:));

tol = 1e-6; % stopping criteria
iter_max = 1000; % maximun number of iterations

hh = h*h;  % h^2

% values of omega
Omega = 1.2:0.01:1.7;

Iter = zeros(size(Omega));  % iteration number for different omega

% start SSOR
for k_omega = 1 : length(Omega)
    omega = Omega(k_omega);
    v = v0;
    for k = 1 : iter_max
        for i = 2:N+1
            for j = 2:N+1
                vij_old = v(i,j);
                v(i,j) = 0.25* (hh*f + v(i+1,j) + v(i-1,j) ...
                    + v(i,j+1) + v(i,j-1));
                v(i,j) = (1-omega)*vij_old + omega*v(i,j);
            end
        end
        for i = N+1:-1:2
            for j = N+1:-1:2
                vij_old = v(i,j);
                v(i,j) = 0.25* (hh*f + v(i+1,j) + v(i-1,j) ...
                    + v(i,j+1) + v(i,j-1));
                v(i,j) = (1-omega)*vij_old + omega*v(i,j);
            end
        end
        
        relerr = norm(v(:) - vt(:)) / norm_vt;
        if relerr < tol
            break;
        end
    end
    Iter(k_omega) = k;
end

% plot the result
plot(Omega,Iter,'o-');
axis([1.2,1.7,20,70])
xlabel('\omega'); ylabel('iteration number');
legend([int2str(N),' x ', int2str(N)]);
title('SSOR with different \omega');


