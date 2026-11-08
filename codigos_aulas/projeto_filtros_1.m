% clear all;
% close all;
% clc;

% % Abre a toolbox de projeto de filtros
% fdatool

fs = 256;                   % Frequência de amostragem
N = 1000;                   % Número de amostras
t = 0:(1/fs):(N-1)*(1/fs);  % Vetor de tempo

y = sin(2*pi*25*t) + sin(2*pi*60*t) + 0.2*randn(1, N); % Sinal de entrada (soma de senóides)

figure(2);
plot(t, y, 'b');
legend('Sinal de entrada');

freqz(SOS, 512, fs);         % Resposta em frequência do filtro
yfilt = filtfilt(SOS, G, y); % Filtragem do sinal de entrada

figure(3);
plot(t, y, 'b', t, yfilt, 'r');
legend('Sinal de entrada', 'Sinal filtrado');