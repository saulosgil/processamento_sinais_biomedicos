% Exemplo 2.9
%
% File neural_data.mat contains two waveforms, x and y, that were recorded from two differ
% ent neurons in the brain with a sampling interval of 0.2 ms. They are believed to be involved in the 
% same function, but separated by one or more neuronal junctions that impart a delay to the signal. 
% Plot the original data, determine if they are related and, if so, the time delay between them.

% Solution

% Take the cross-correlation between the two signals using axcor. Find the maximum correla
% tion and the time shift at which that maximum occurs. The former will tell us if they are related 
% and the latter the time delay between the two nerve signals.
clear all;
close all;
clc;

addpath(genpath("funcoes_auxiliares_and_data_tests"));   % torna axcor visível

caminho = "funcoes_auxiliares_and_data_tests\Chapter 2\neural_data.mat";
load(caminho);

fs = 1/0.0002;                        % frequência de amostragem (1/Ts)
t = (1:length(x))/fs;                 % vetor de tempo

subplot(2,1,1);
plot(t, y, 'k', t, x, ':');           % dados originais
xlabel('Tempo (s)');
ylabel('Amplitude');
legend('y(t)', 'x(t)');

[rxy, lags] = axcor(x, y);            % correlação cruzada
% alternativa sem a função do livro:
% [rxy, lags] = xcorr(x - mean(x), y - mean(y), 'coeff');

[max_corr, idx] = max(rxy);           % correlação máxima e seu índice
max_shift = lags(idx)/fs;             % deslocamento em segundos

subplot(2,1,2);
plot(lags/fs, rxy, 'k'); hold on;     % função de correlação cruzada
plot(max_shift, max_corr, '*r', 'MarkerSize', 10);   % ponto de máximo
hold off;
xlabel('Lag (s)');
ylabel('r_{xy}');
legend('Correlação cruzada', 'Máximo');

fprintf('Correlação máxima: %.4f  |  Deslocamento: %.4f s\n', max_corr, max_shift);

% Análise
%
% Após a correlação cruzada, encontrar a correlação máxima é simples usando o 
% operador max do MATLAB. 
% Encontrar o instante em que esse valor máximo ocorre é 
% um pouco mais complicado. O operador max fornece o índice do valor máximo, 
% chamado aqui de max_shift. Para encontrar o deslocamento real correspondente a 
% esse índice, é preciso obter o valor de lag nessa posição, isto é, 
% lags(max_shift). 
% Esse valor de lag então precisa ser convertido no deslocamento 
% temporal correspondente, dividindo-o pela frequência de amostragem, fs (ou 
% multiplicando-o por Ts). A posição do pico é marcada na curva de correlação 
% cruzada como um ponto * na Figura 2.14.

% Resultado
%
% Os dois sinais são mostrados na Figura 2.14a e a função de correlação cruzada aparece
% na Figura 2.14b (ou ao rodar esse script). O instante em que o pico ocorre está indicado no
% gráfico e, após a conversão para tempo descrita acima, corresponde a 0,013 s. 
% A correlação máxima é 0,45, o que sugere que os dois sinais estão relacionados.

% O próximo exemplo antecipa a transformada de Fourier, apresentada no Capítulo 3. 
% No Exemplo 2.10, pegamos o sinal de EEG mostrado na Figura 2.15 e o comparamos com 
% uma senoide de determinada frequência. Como o MATLAB fará o trabalho, podemos realizar 
% essa comparação ao longo de uma faixa de frequências. Infelizmente, não podemos simplesmente
% comparar o sinal com uma onda seno usando a correlação padrão (por exemplo, a Equação 2.29), 
% porque queremos compará-lo com uma senoide genérica, e não apenas com uma onda seno. 
% Uma senoide pode ser uma onda seno, mas com qualquer deslocamento de fase (ou deslocamento temporal
% equivalente) de até 90°, inclusive. 
% Uma abordagem é usar a correlação cruzada para comparar o sinal de EEG com ondas seno deslocadas e
%  tomar a correlação máxima, como no exemplo anterior. (Se fôssemos espertos, limitaríamos o deslocamento
% ao equivalente a ±90°, mas é mais fácil usar o deslocamento máximo e deixar o MATLAB fazer o trabalho 
% extra.)