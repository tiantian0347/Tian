function [x, relres, iter, flag, resvec] = PCG(A, x, b, M, IterMax, tol)

% 预处理共轭梯度法
% 输入参数:
%   A: 系数矩阵
%   x: 迭代初值
%   b: 右端项
%   M: 预处理子
%   IterMax: 最大迭代步数
%   tol: 停机准则
%
% 输出参数:
%   x: 近似解
%   relres: 相对残量范数
%   iter: 迭代步数
%   flag: 算法成功标识
%   resvec: 每次迭代的相对残量
%


  flag = 0;
  iter = 0;

  alpha = 0.0;
  beta  = 0.0;
  bnrm2 = 0.0;
  error = 0.0;
  rho   = 0.0;
  rho_1 = 0.0;

  [dim,dim] = size(A);

  p     = zeros(dim,1);
  q     = zeros(dim,1);
  r     = zeros(dim,1);
  z     = zeros(dim,1);

  % -----------------------------
  % Quick check of approximation.
  % -----------------------------

  bnrm2 = norm( b );
  if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end

  r = b - A*x;
  error = norm( r ) / bnrm2;
  if ( error < tol ) return, end

  % ----------------
  % Begin iteration.
  % ----------------

  for iter = 1:max_it

     z  = M \ r;
     rho = (r'*z);

     % -------------------------
     % Compute direction vector.
     % -------------------------

     if ( iter > 1 ),
        beta = rho / rho_1;
        p = z + beta*p;
     else
        p = z;
     end

     q = A*p;
     alpha = rho / (p'*q );

     % ---------------------
     % Update approximation.
     % ---------------------

     x = x + alpha * p;

     r = r - alpha*q;

     % ------------------
     % Check convergence.
     % ------------------

     error = norm( r ) / bnrm2;
     if ( error <= tol ), break, end 

     rho_1 = rho;

  end

  % ------------------------
  % Final convergence check.
  % ------------------------

  if ( error > tol ) flag = 1; end

% --------
% End cg.m
% --------
