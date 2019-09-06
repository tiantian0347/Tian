%% Least squares smoothing 
% This example illustrates smooth a noisy ECG waveform
% by least squares.
% This approach regularizes the energy of the second-order derivative.
%
%  Ivan Selesnick
% selesi@poly.edu

%% Start

clear
close all

%% Load data

load data.txt;
whos

y = data;         % data value
N = length(y);

%% Display data

figure(1)
clf
plot(y)
title('Data')

%% Smoothing (degree = 2)
% D is the second-order difference matrix.
% It approximates the second-order derivative.
% In order to exploit fast banded solvers in Matlab,
% we define D as a sparse matrix using 'spdiags'.

e = ones(N, 1);
D = spdiags([e -2*e e], 0:2, N-2, N);

%%
% Observe the first and last corners of D.
% (D is too big to display, so we show
% the first and last corners only.)

%%
% First corner of D:

full(D(1:5, 1:5))

%%
% Last corner of D:

full(D(end-4:end, end-4:end))

%%
% Solve the least square problem

lam = 50;
F = speye(N) + lam * D' * D;            % F is a banded matrix
x = F \ y;                              % Matlab uses a fast solver for banded systems)

figure(1)
plot(x)

