% Least squares smoothing 
% This example illustrates smooth a noisy ECG waveform by least squares.
% This approach regularizes the energy of the second-order derivative.
%
%  Ivan Selesnick
% selesi@poly.edu

clear all;
close all;

load data.txt;
y = data;        
N = length(y);

figure(100)
plot(y)
title('observed data')

% Smoothing (degree = 2)
% D is the second-order difference matrix.
% It approximates the second-order derivative.
% In order to exploit fast banded solvers in Matlab,
% we define D as a sparse matrix using 'spdiags'.

e = ones(N, 1);
D = spdiags([e -2*e e], 0:2, N-2, N);

% Observe the first and last corners of D
full(D(1:5, 1:5))
full(D(end-4:end, end-4:end))

% Solve the least square problem
Alpha = [0.1,0.5,1,5,10,50];
for k = 1 : length(Alpha)
    alpha = Alpha(k);
    F = speye(N) + alpha * D' * D;
    x = F \ y;      % Matlab uses a fast solver for banded systems
    
    figure(k); plot(x); title(['\alpha=',num2str(alpha)]);
    pause(1)
end

