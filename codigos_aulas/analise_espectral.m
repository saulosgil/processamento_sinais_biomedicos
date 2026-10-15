% Verifica análise espectral - Influência das Janelas
% 19/10/2026
% Saulo Gil
clear all;
close all;
clc;

addpath(genpath("funcoes_auxiliares_and_data_tests"));   % torna axcor visível

% Gerando sinais
fs = 1000;
N = 1024;
x = sig_noise([235 250], -3, N);
f = (0:N-1)*(fs/N);
t = (0:N-1)*(1/fs);

figure(1);
plot(t, x);
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal no tempo')
grid on;

Xmag = abs(fft(x)); % espectro de magnitude
figure(2);
subplot(3,1,1);
plot(f(1:N/2), Xmag(1:N/2));
xlabel('Freq (Hz)');
ylabel('|X(f)|');
title('Espectro de Magnitude')
grid on;

% Janela de Hamming
x1 = x.*hamming(N)';
Xmag1 = abs(fft(x1));
figure(2); hold on; 
subplot(3,1,2);
plot(f(1:N/2), Xmag1(1:N/2));
xlabel('Freq (Hz)');
ylabel('|X(f)|');
title('Espectro de Magnitude - Janela de Hamming')
grid on;

% Janela de BlackmanHarris
x2 = x.*blackmanharris(N)';
Xmag2 = abs(fft(x2));
figure(2); hold on; 
subplot(3,1,3);
plot(f(1:N/2), Xmag2(1:N/2));
xlabel('Freq (Hz)');
ylabel('|X(f)|');
grid on;
title('Espectro de Magnitude - Janela de Blackman-Harris')
