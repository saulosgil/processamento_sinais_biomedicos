% Manipulação de dados e introdução ao MATLAB
% 01/10/2026
% Saulo Gil

clear all; % limpa o wordspace
close all; % fecha todas as figuras 

load("funcoes_auxiliares_and_data_tests\Chapter 2\eeg_data.mat"); % carrega os dados do EEG

figure(1)
plot(eeg); % plota o primeiro canal do EEG
xlabel('amostras');
ylabel('amplitude');
title('Dados de EEG');
xlim([0 801]); % limita o eixo x para as primeiras 801 amostras
grid on; % adiciona uma grade ao gráfico

% calculo da média

Namostras = length(eeg); % calcula o número de amostras
media = (1 / Namostras) * sum(eeg); % calcula a média do sinal
media

media_fct = mean(eeg); % calcula a média do sinal
media_fct

% calculo da variancia e desvio padrão
variancia = (1 / (Namostras - 1)) * sum((eeg - media).^2); % calcula a variância do sinal
variancia

variancia_fct = var(eeg); % calcula a variância do sinal
variancia_fct

desvio_padrao = sqrt(variancia); % calcula o desvio padrão do sinal
desvio_padrao

desvio_padrao_fct = std(eeg); % calcula o desvio padrão do sinal
desvio_padrao_fct

% Lembrete-----------------------
% usar help nome_da_funcao para obter informações sobre a função
% ex. help mean
% help mean 

% Calculo do RooT Mean Square (RMS)
rms = sqrt(mean(eeg.^2)); % calcula o RMS do sinal
rms

% Calculo da energia do sinal
energia = sum(eeg.^2); % calcula a energia do sinal         
energia

% Calculo da potência do sinal
potencia = energia / Namostras; % calcula a potência do sinal
potencia

% Normalização do sinal - Z-Score
eeg_normalizado = (eeg - media) / desvio_padrao; % normaliza o sinal
eeg_normalizado

media_normalizada = mean(eeg_normalizado); % calcula a média do sinal normalizado
media_normalizada

desvio_padrao_normalizado = std(eeg_normalizado); % calcula o desvio padrão do sinal normalizado
desvio_padrao_normalizado

% Plot do sinal normalizado
figure(2)
plot(eeg_normalizado); % plota o sinal normalizado
xlabel('amostras');
ylabel('amplitude');
title('Sinal Normalizado (Z-Score)');
xlim([0 801]); % limita o eixo x para as primeiras 801 amostras
grid on;

% Características do sinal normalizado
figure(3)
data_histogram = histogram(eeg_normalizado); % plota o histograma do sinal normalizado
grid on;

figure(4)
data_histogram_densidadeProb = histogram(eeg_normalizado, 'Normalization', 'pdf'); % plota o histograma do sinal normalizado com normalização de densidade de probabilidade
grid on;

% Função de densidade de probabilidade (PDF) do sinal normalizado
y = -5:0.05:5; % cria um vetor de valores para o eixo x
mu = mean(eeg_normalizado); % calcula a média do sinal normalizado
sigma = std(eeg_normalizado); % calcula o desvio padrão do sinal normalizado

f = exp(- (y - mu).^2./ (2 * sigma^2)) ./ (sigma * sqrt(2 * pi)); % calcula a função de densidade de probabilidade do sinal normalizado

figure(4); hold on; % mantém o gráfico atual
plot(y, f, 'LineWidth', 1.5); % plota a função de densidade de probabilidade do sinal normalizado
hold off; % libera o gráfico atual


% Exemplo 1.5
% Senoide de 500 pontos com frequência de amostragem de 500 Hz e frequencia de 4 Hz.
% Calcule
fs = 500; % frequência de amostragem
N = 500; % número de pontos
f = 4; % frequência do sinal

t = (0:N-1)/fs; % vetor de tempo
sinal = sin(2*pi*f*t); % gera o sinal senoide 

figure(5)
plot(t, sinal); % plota o sinal senoide

rms_sinal = sqrt(mean(sinal.^2)); % calcula o RMS do sinal senoide
rms_sinal  % valor conhecido do RMS de uma senoide é 0.7071, como o do calculado acima.



