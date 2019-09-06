% 二维离散 Poisson 方程的 ILU 分解
%
clear all;
N = 4;
T1D = 2*eye(N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1);
T2D = kron(eye(N),T1D) + kron(T1D,eye(N));
T2D = sparse(T2D);

[L,U] = lu(T2D);
setup.type = 'nofill';
[IL,IU] =  ilu(T2D,setup);

figure(11); spy(T2D);
figure(12)
subplot(2,2,1); spy(L);
subplot(2,2,2); spy(U);
subplot(2,2,3); spy(IL);
subplot(2,2,4); spy(IU);

fprintf('||T2D-L*U||_F=%.4e\n', norm(T2D-L*U,'fro'))
fprintf('||T2D-IL*IU||_F=%.4e\n', norm(T2D-IL*IU,'fro'))

