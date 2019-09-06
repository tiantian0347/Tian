function [beta,v] = House(x)
%
% Householder transformation
%
% Given x, compute v such that Hx=\|x\|_2 e_1 where H = I - beta*v*v^T
%

n = length(x);
sigma = dot(x(2:n),x(2:n));
v = x;
if sigma == 0 
    if x(1) < 0
        v(1) = 2*x(1);
        beta = 2/(v(1)*v(1));
    else
        v(1) = 1; 
        beta = 0;
    end
else
    alpha = sqrt(x(1)*x(1)+sigma);
    if x(1) < 0
        v(1) = x(1) - alpha;
    else
        v(1) = -sigma/(x(1)+alpha);
    end
    beta = 2/(v(1)*v(1)+sigma);
end
