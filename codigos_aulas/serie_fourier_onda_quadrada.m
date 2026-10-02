% Implementa Série de Fourier para uma onda quadrada.
% 02/10/2026
% Saulo Gil
% A série de Fourier é uma representação de funções periódicas como uma soma infinita de senos e cossenos. 
% Para uma onda quadrada, a série de Fourier pode ser expressa como:
% f(t) = (4/π) * Σ (1/n) * sin(nωt), onde n é um número ímpar (1, 3, 5, ...), 
%     ω é a frequência angular e 
%     t é o tempo.

close all;
clear all;

t = -4:0.01:4; % Intervalo de tempo
y = zeros(length(t), 1); % Inicializa o vetor de saída

% Define a onda quadrada
for kk = 1:length(t)
    if t(kk) >= -pi/2 & t(kk) <= pi/2
        y(kk) = 1;
    end
end 

% Plot da onda quadrada
figure(1);
plot(t, y, 'r', 'LineWidth', 2);
hold on;
grid on;
xlabel('Tempo');
ylabel('Amplitude');
title('Intro á Série de Fourier - Onda Quadrada');
axis([-4 4 -0.5 1.5]);

% Criando aproximação da Série de Fourier
% começando pela média da função

x = (1/2)*ones(1, length(t)); % Inicializa a aproximação da série de Fourier

% Loop para mostrar a aproximação da série de Fourier com diferentes números de termos
for n = 1:2:11
    x = x + (2/(n*pi))*sin(n*pi/2)*cos(n*t); % Adiciona os termos da série de Fourier
    figure(1);
    hold on;
    plot(t, x, 'b'); % Plota a aproximação da série de Fourier
end

% -- Comentário --
% Repare que a aproximação da série de Fourier melhora à medida que mais termos são adicionados.

% Loop para mostrar a aproximação da série de Fourier com diferentes números de termos isoladamente
% começando pela média da função
% Descomentar a linha abaixo para ver a aproximação com 101 termos
% x = (1/2)*ones(1, length(t)); % Inicializa a aproximação da série de Fourier

% for n = 1:2:101
%     x = x + (2/(n*pi))*sin(n*pi/2)*cos(n*t); % Adiciona os termos da série de Fourier
% end

% figure(1);
% hold on;
% plot(t, x, 'b'); % Plota a aproximação da série de Fourier

% -- Comentário --
% Repare que a aproximação da série de Fourier melhora à medida que mais termos são adicionados.
% Aqui foram adicionados 101 termos, o que resulta em uma aproximação muito próxima da onda quadrada original.








