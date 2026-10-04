% Exemplo 2.7
% Find the Pearson correlation coefficient between the two waveforms shown in Figure 2.13. 
% These waveforms are stored as variables x and y in file Ex2_7.mat.

% Solution
% Load the file and find the unnormalized correlation using Equation 2.29 as in Example 2.5 and 
% apply Equation 2.44. Subtract the means of each signal before correlation.
clear all;
close all;
clc;

caminho = "funcoes_auxiliares_and_data_tests\Chapter 2\Ex2_7.mat"; % caminho dos dados
load(caminho);                                                     % carrega os dados do arquivo .mat

N = length(x);                             % Encontra N
rxy = sum((x - mean(x)) .* (y - mean(y))); % subtrai a média e aplica a Eq. 2.29
rxy = rxy / ((N-1) *sqrt(var(x)*var(y)));  % Aplica a Eq. 2.44
disp(['Correlation: ', num2str(rxy)]);     % Output correlation

% Comparação com funções nativas
Rxx = corrcoef(x, y); % calcula a matrix de correlação 
disp(['Correlation com função nativa: ', num2str(Rxx(1,2))]);    

% calculo da matrix de covariância
S = cov(x, y); % calcula a matrix de covariância
S