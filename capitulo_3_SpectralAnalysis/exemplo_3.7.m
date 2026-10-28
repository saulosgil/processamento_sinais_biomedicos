% Example 3.7
% Evaluate the influence of averaging power spectra on a signal made up of a combination of 
% broadband and narrowband signals along with added noise. The data can be found in file 
% broadband1.mat, which contains a white-noise signal filtered to be between 0 and 300 Hz 
% and two closely spaced sinusoids at 390 and 410 Hz.

% Solution
% Load the test data file broadband1 containing the narrowband and broadband processes. 
% First, calculate and display the unaveraged PS using Equation 3.34. Then apply PS averaging 
% using an averaging routine. Use a segment length of 128 points with a maximum overlap of 127 
% points. To implement averaging, use pwelch with the appropriate calling parameters.
clear all;
close all;
clc;

addpath(genpath("funcoes_auxiliares_and_data_tests"));   % torna broadband1.mat visível

% Example 3.7 Investigation of the use of averaging to improve broadband spectral characteristics in the power spectrum.
load broadband1;                            % load data 

fs = 1000;                                  % sampling frequency
nfft = 128;                                 % segment lenght

% Un-averaged spectrum using the direct method of Eq. 3.34

PS = abs((fft(x)).^2)/length(x);            % calculate un-average PS
half_length = fix(length(PS)/2);            % valid points
f = (0:half_length-1)*fs/(2*half_length);   % frequency vector for plotting

subplot(1,2,1)
plot(f,PS(1:half_length),'k'); 

[PS_avg,f] = pwelch(x,nfft, nfft-1,[], fs); % Periodogram, max. overlap

subplot(1,2,2)
plot(f,PS_avg,'k'); 


% Não foi possível replicar pois o arquivo contento os dados (broadband1) não está disponível
