% Example 2.12 Convolution of a first-order system with a step
%
% Convolução é a operação matemática que calcula a saída de um sistema linear e invariante no tempo a 
% partir da entrada e da resposta ao impulso. Se você conhece h(t), consegue prever o que o sistema faz 
% com qualquer sinal.

% A ideia por trás

% Todo sinal pode ser visto como uma sequência de impulsos, um em cada instante, cada um com a amplitude 
% do sinal naquele ponto. O sistema responde a cada impulso com uma cópia de h(t), escalada pela amplitude 
% daquele impulso e deslocada para o instante em que ele ocorreu. Como o sistema é linear, a saída é a soma 
% de todas essas respostas sobrepostas.

% No exemplo da exponencial, como na figura deste exemplo, cada amostra da entrada “dispara” uma pequena 
% exponencial decrescente. A saída em um instante qualquer é a soma das caudas de todas as exponenciais 
% disparadas antes dele. Por isso o sistema tem memória: o que entrou há pouco ainda pesa bastante, e o 
% que entrou há muito tempo quase não pesa.

clear all;
close all;
clc;

fs = 500; 
N = 2500; 
t = (0:N-1)/fs; 
tau = 1; 
h = exp(-t./tau); 
x = ones(1,N); 
y = conv(x,h);

subplot(1,2,1);
plot(t,h); 
title('Impulse Response');
xlabel('Time (s)');
ylabel('h(t)');
grid on;

subplot(1,2,2);
plot(t,y(1:N));
title('Step Response');
xlabel('Time (s)');
ylabel('y(t)');
grid on;

% Resultados
%
% As respostas ao impulso e ao degrau de um sistema de primeira ordem são mostradas na Figura. 
% Observe que a função conv do MATLAB gera várias amostras adicionais (na verdade, 4999 pontos); 
% portanto, apenas os primeiros 2500 pontos são utilizados na representação gráfica da resposta.

% No processamento de sinais, a convolução pode ser utilizada para implementar os filtros básicos 
% descritos no Capítulo 4. Assim como seus equivalentes analógicos, os filtros digitais são processos 
% lineares que modificam o espectro do sinal de entrada de alguma maneira desejada, por exemplo, para 
% reduzir o ruído. Como ocorre com todos os processos lineares, a resposta ao impulso do filtro, 
% h[k] ou h(t), descreve completamente o filtro.

% Leitura prática

% Juntos, os dois gráficos mostram um sistema que suaviza a entrada: mudanças bruscas aparecem na saída 
% de forma gradual, com atraso controlado por τ. Isso equivale a um filtro passa-baixas, com frequência 
% de corte fc = 1/(2πτ) ≈ 0,16 Hz neste caso. Em biossinais, esse comportamento aparece, por exemplo, na 
% resposta de eletrodos, em filtros RC dos amplificadores e na dinâmica de vários processos fisiológicos 
% lentos.