% Evaluate a signal to determine if it is stationary. 
% If not, attempt to modify the signal to remove the nonstationarity, if possible. 
% The file data_c1.mat contains the signal in variable x. 
% The signal is 1000 points long and the sample interval is Ts = 0.001 s.

clear all; 
close all; 
clc;

% Example 2.2
% Part 1. Evaluate a signal for stationarity.
caminho = "funcoes_auxiliares_and_data_tests\Chapter 2\data_c1.mat"; % carrega os dados
load(caminho); % carrega os dados do arquivo .mat

for k = 1:4                       % Loop através de quatro segmentos de 250 pontos cada
    m = 250*(k-1) + 1;            % Define o índice inicial do segmento  
    segment = x(m:m + 249);       % Define o segmento de 250 pontos do sinal x  
    avg(k) = mean(segment);       % Calcula a média do segmento
    variance(k) = var(segment);   % Calcula a variância do segmento
end

disp('Mean Segment 1 Segment 2 Segment 3 Segment 4')      % Heading
disp(avg)                                                 % Apresenta a média de cada segmento

disp('Variance Segment 1 Segment 2 Segment 3 Segment 4')  % Heading
disp(variance)                                            % Apresenta a variância de cada segmento
disp('The signal is nonstationary because the mean and variance are not constant.') % Mensagem indicando que o sinal é não estacionário
disp('------------------------------------------------------------------------------')

%%
% Example 2.8
% Part 2. Modify a nonstationary signal to become stationary
%
y = [diff(x);0]; % calcula a primeira diferença do sinal x e adiciona um zero no final para manter o mesmo comprimento

for k = 1:4                       % Loop através de quatro segmentos de 250 pontos cada
    m = 250*(k-1) + 1;            % Define o índice inicial do segmento  
    segment = y(m:m + 249);       % Define o segmento de 250 pontos do sinal y  
    avg(k) = mean(segment);       % Calcula a média do segmento
    variance(k) = var(segment);   % Calcula a variância do segmento
end
 
disp('Mean Segment 1 Segment 2 Segment 3 Segment 4')      % Heading
disp(avg)                                                 % Apresenta a média de cada segmento

disp('Variance Segment 1 Segment 2 Segment 3 Segment 4')  % Heading
disp(variance)                                            % Apresenta a variância de cada segmento
disp('The signal is now stationary because the mean and variance are constant.') % Mensagem indicando que o sinal é estacionário
disp('------------------------------------------------------------------------------')