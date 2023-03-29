function [Phi, omega, lambda, b, freq, Xdmd ,r] = DynamicModeDecomp(X, varin)
% function [Phi,omega,lambda,b, freq, Xdmd] = DynamicModeDecomp(X,varin)
% Computes the Dynamic Mode Decomposition of data X
%
% INPUTS: 
% Columns of X are state snapshots 
% Rows of X are measurements
% Columns X are time points, sampled at equal dt's
%
% Optional parameters: {'parameter_name', [default_value]}
%   {'dt', [1]}
%   {'r', [1e32]}        truncate to rank-r
%   {'nstacks', 1}       number of stacks of the raw data
%
%
% OUTPUTS:
% Phi, the DMD modes
% omega, the continuous-time DMD eigenvalues
% lambda, the discrete-time DMD eigenvalues
% b, a vector of magnitudes of modes Phi
% freq, the estimated frequencies in Hertz
% Xdmd, the data matrix reconstrcted by Phi, omega, b
%
%
% References:
% Kutz, J. N., Brunton, S. L., Brunton, B. W., & Proctor, J. L. (2016). 
% Dynamic Mode Decomposition. Society for Industrial and Applied 
% Mathematics. 
% https://doi.org/10.1137/1.9781611974508
%
% Note: This code is based on DMDfull.m provided by these authors. The only
% changes are the addition of frequency output and the time range change in
% the calculation of Xdmd. 
%
%
% Efe Ilicak, 30/10/2022.



%% input parsing
p = inputParser; 

% required inputs
p.addRequired('X', @isnumeric);

% parameter value iputs
p.addParameter('dt', 1, @(x)isnumeric(x) && x>0);
% p.addParameter('r', 1e32, @(x)isnumeric(x) && x>0);
p.addParameter('nstacks', 1, @(x)isnumeric(x) && x>0);

% now parse the inputs
p.parse(X, varin{:});
inputs = p.Results;

%% stacking the data matrix 
if inputs.nstacks > 1
    Xaug = [];
    for st = 1:inputs.nstacks
        Xaug = [Xaug; X(:, st:end-inputs.nstacks+st)];
    end
    X1 = Xaug(:, 1:end-1);
    X2 = Xaug(:, 2:end);
else
    X1 = X(:, 1:end-1);
    X2 = X(:, 2:end);
end
m = size(X1, 2);

%% DMD
[U, S, V] = svd(X1, 'econ');

% b = diff(log(diag(S)));
% [~,loc] = min(b(2:end));
% r = loc+2;
r = 15;

% r = inputs.r;
if r >= size(U, 2)
    r = size(U, 2);
end

if r >= size(U,2) % no rank truncation
    Atilde = U' * X2 * V / S;
    [W, D] = eig(Atilde); 
    Phi = X2 * V / S * W;
else % truncate modes
    U_r = U(:, 1:r);
    S_r = S(1:r, 1:r);
    V_r = V(:, 1:r);
    Atilde = U_r' * X2 * V_r / S_r;
    [W_r, D] = eig(Atilde);
    Phi = X2 * V_r / S_r * W_r;
end

lambda = diag(D);
omega = log(lambda)/inputs.dt;

%% compute DMD mode amplitudes
x1 = X1(:, 1);
b = Phi\x1;

%% compute the frequencies
freq = atan2(imag(lambda),real(lambda))/(2*pi*inputs.dt);

%% DMD reconstructions
time_dynamics = zeros(r, m);
t = (0:m)*inputs.dt;
for iter = 1:m
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
Xdmd = Phi * time_dynamics;

