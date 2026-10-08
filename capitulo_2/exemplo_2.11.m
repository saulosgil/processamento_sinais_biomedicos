% Example 2.11 Use of autocovariance to determine the correlation
% of heart rate variation between heart beats
%
clear all;
close all;
clc;

addpath(genpath("funcoes_auxiliares_and_data_tests"));   % torna axcor visível

caminho = 'funcoes_auxiliares_and_data_tests\Chapter 2\Hr_pre.mat';
load(caminho);

figure(1)
plot(hr_pre) % plot para ver o sinal

[cov_pre,lags_pre] = axcor(hr_pre - mean(hr_pre));  % Auto-covariance
figure(2)
plot(lags_pre,cov_pre,'k');                         % Plot normal auto-cov
grid on;
hold on; 

plot([lags_pre(1) lags_pre(end)], [0 0],'k');       % Plot a zero line
axis([-30 30 -0.2 1.2]);                            % Limit x-axis to ± 30 beats

% Results
%
% Os resultados da Figura mostram que existe correlação entre batimentos cardíacos 
% adjacentes ao longo de até 10 batimentos. A autocovariância dos dados da frequência 
% cardíaca durante a meditação é examinada no Problema 2.35.


