% Executa visualização de espectro de frequencia de uma dada função de transferência de filtro digital.
close all;
clear all;
clc;

fs = 1000;  % Frequência de amostragem
N = 512;    % Número de pontos para a FFT

a = [1 -0.2 0.8];  % Coeficientes do denominador (filtro FIR)
b = [0.2 0.5]; % Coeficientes do numerador (filtro FIR)

% Resposta em frequência via DFT dos coeficientes do filtro
H = fft(b,N)./fft(a,N);           % espectro de frequencia do filtro de H
Hmag = 20*log10(abs(H));          % Magnitude do espectro
Theta = angle(H)*360/(2*pi);      % Fase do espectro em graus

f = (1:N)*fs/N;                   % Eixo de frequencia

% Plota a magnitude do espectro de frequencia
figure;
subplot(2,1,1);
plot(f(1:N/2), Hmag(1:N/2),'LineWidth',2);
grid on;
title('Magnitude do Espectro de Frequência');
xlabel('Frequência (Hz)');
ylabel('Magnitude (dB)');

subplot(2,1,2);
plot(f(1:N/2), Theta(1:N/2),'LineWidth',2);
grid on;
title('Fase do Espectro de Frequência');
xlabel('Frequência (Hz)');
ylabel('Fase (graus)');

% Analisando resposta impusiva
x = [1, zeros(1, N-1)];     % Sinal de entrada (impulso unitário)
y = filter(b,a,x);          % Resposta ao impulso
t = 0:(1/fs):(N-1)*(1/fs);  % Eixo de tempo

figure;
plot(t, y, 'LineWidth', 2);
grid on;
title('Resposta ao Impulso');
xlabel('Tempo (s)');
xlim([0 0.1]);
ylabel('Magnitude');    

% Analisando resposta em frequência usando freqz
% Passa os coeficientes do filtro para a função freqz, que calcula e plota a resposta em frequência do 
% filtro digital.
freqz(b, a, N, fs);