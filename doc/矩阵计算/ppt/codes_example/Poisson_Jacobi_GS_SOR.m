%
% Jacobi, GS and SOR for the 2-D Poisson problem
% discretized by five point finite difference scheme
%

clear all;
close all;
N = 32;
h = 1/(N+1); % step size
v0 = zeros(N+2,N+2); % initial guess + boundary condition
x = 0:h:1;  y = 0:h:1;
v0(1,:) = 0.25 * (y .* y);
v0(N+2,:) = 0.25 * (1 + y .* y);
v0(:,1) = 0.25 * (x .* x);
v0(:,N+2) = 0.25 * (x .* x + 1);

f = -1; % f(x,y)
[X,Y] = meshgrid(x);
vt = 0.25*(X.*X + Y.*Y); % true solution
norm_vt = norm(vt(:));

iter = 200; % number of iteration steps

hh = h*h;

% Jacobi iteration
relerr_Jacobi = zeros(iter+1,1); % relative error
relerr_Jacobi(1) = 1;

v = v0;
for k = 1 : iter
    v_old = v;
    for i = 2:N+1
        for j = 2:N+1
            v(i,j) = 0.25* (hh*f + v_old(i+1,j) + v_old(i-1,j) ...
                + v_old(i,j+1) + v_old(i,j-1));
        end
    end
    relerr_Jacobi(k+1) = norm(v(:) - vt(:)) / norm_vt;
end

% GS iteration
relerr_GS = zeros(iter+1,1); % relative error
relerr_GS(1) = 1;

v = v0;
for k = 1 : iter
    for i = 2:N+1
        for j = 2:N+1
            v(i,j) = 0.25* (hh*f + v(i+1,j) + v(i-1,j) ...
                                 + v(i,j+1) + v(i,j-1));
        end
    end
    relerr_GS(k+1) = norm(v(:) - vt(:)) / norm_vt;
end

% SOR iteration
relerr_SOR = zeros(iter+1,1); % relative error
relerr_SOR(1) = 1;

omega = 2 / (1+sin(pi/(N+1)));
v = v0;
for k = 1 : iter
    for i = 2:N+1
        for j = 2:N+1
            vij_old = v(i,j);
            v(i,j) = 0.25* (hh*f + v(i+1,j) + v(i-1,j) ...
                                 + v(i,j+1) + v(i,j-1));
            v(i,j) = (1-omega)*vij_old + omega*v(i,j);
        end
    end
    relerr_SOR(k+1) = norm(v(:) - vt(:)) / norm_vt;
end

% plot the error
step = 4;
xx = 1 : step : iter+1;
semilogy(xx, relerr_Jacobi(xx),'-b+');
hold on
semilogy(xx, relerr_GS(xx),'-ro');
semilogy(xx, relerr_SOR(xx),'-kd');
legend('Jacobi', 'GS', 'SOR');
title(['N=',int2str(N)]);

