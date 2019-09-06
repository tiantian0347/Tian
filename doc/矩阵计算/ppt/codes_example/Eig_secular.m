% 对称矩阵特征值 分而治之算法 特征方程

clear all; 
close all;

alpha=0.005;
f = @(x) 1./(1-x) + 1./(2-x) + 1./(3-x) + 1./(4-x);

x = -1:0.0001:6;
y = 1 + alpha*f(x);
plot(x,y,'LineWidth',1)
axis([-1,6,-6,6]);

% hold on
% plot(x,ones(size(x)),'g-.')
% plot(x,zeros(size(x)),'g-.')
% t = -6:0.001:6;
% plot(ones(size(t)), t, 'g-.')
% plot(2*ones(size(t)), t, 'g-.')
% plot(3*ones(size(t)), t, 'g-.')
% plot(4*ones(size(t)), t, 'g-.')

tit = ['\alpha=',num2str(alpha)];
title(tit,'FontSize',16)
set(gca,'FontSize',16)
% print('-dpng','secular.png')
