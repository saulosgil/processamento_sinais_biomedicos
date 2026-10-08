% Exemplo 2.8

% Determine if a sine wave and cosine wave at the same frequency are orthogonal and if sine waves 
% at harmonically related frequencies are orthogonal. The term “harmonically related” means 
% that sinusoids are related by frequencies that are multiples. Thus, the signals sin(2t), sin(4t), 
% and sin(6t) are harmonically related.

% Solution

% Generate a 500-point, 1.0-s time vector as in the last example. Use this time vector to generate a 
% data matrix where the columns represent 2- and 4-Hz cosine and sine waves. 
% Apply the covariance and correlation MATLAB routines (i.e., cov and corrcoef) and display results.

% Example 2.8 Application of the covariance matrix to sinusoids that are orthogonal and a sawtooth
clear all;
close all;
clc;
%
N = 1000;              % Number of points
Tt = 2;                % desired total time
fs = N/Tt;             % Calculate sampling frequency
t = (0:N-1)/fs;        % Time vector
X(:,1) = cos(2*pi*t)'; % Generate a 1 Hz cosine
X(:,2) = sin(2*pi*t)'; % Generate a 1 Hz sine
X(:,3) = cos(4*pi*t)'; % Generate a 2 Hz cosine
X(:,4) = sin(4*pi*t)'; % Generate a 1 Hz sine
%

S = cov(X)             % Print covariance matrix
Rxx = corrcoef(X)      % and correlation matrix

% Resultado
%
% The off-diagonals of both matrices are zero, showing that all these sine and cosine waves are 
% uncorrelated and hence orthogonal.