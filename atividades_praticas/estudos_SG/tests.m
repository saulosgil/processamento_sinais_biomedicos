% Atividade filtragem digital de sinais
% Objetivo:
% Tentar isolar alguma baixa frequência após inspeção da PSD
clear all;
close all;
clc;

addpath(genpath("atividades_praticas"));   % torna broadband1.mat visível

load('S03.mat'); % Carrega os filtros projetados
load('LowPassFilter.mat'); % Carrega os filtros projetados
load('BandPassFilter.mat'); % Carrega os filtros projetados

[Nsamples, Nchannels,trials] = size(stimuliData); % Obtém o número de amostras e canais do sinal
t = 0:(1/fs):(Nsamples-1)*(1/fs);  % Vetor de tempo
sinal = stimuliData(:, 3, 3) % POz, Trial 3, estimulo 10Hz

figure(1);
plot(t, sinal)
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal original no canal POz - Trial 3');
grid on;
 
[pxx, fw] = pwelch(sinal, 1536, fs);

% PSD do sinal original
figure(2);
plot(fw, pxx, 'k');
xlim([1.5 2.5]);
xlabel('Frequência (Hz)');
ylabel('PSD');
title('Densidade espectral de potência do sinal original');
grid on;

% Aplicar o filtro projetado ao sinal
% LowPass - IRR ButterWorth
% Specifity order = 6
% Fs = 256
% Fc = 1.7
sinalFiltrado = filtfilt(SOS, G, sinal);

[pxxf, fwf] = pwelch(sinalFiltrado, 1536, fs);

% PSD do sinal filtrado
figure(3);
plot(fwf, pxxf, 'k');
xlim([0 2.5]);
xlabel('Frequência (Hz)');
ylabel('PSD');
title('Densidade espectral de potência do sinal filtrado - LowPass');
grid on;

% Sinal Filtrado 
figure(4);
plot(t, sinalFiltrado)
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal filtrado (LowPass) no canal POz - Trial 3');
grid on;

% Aplicar o filtro projetado ao sinal
% BandPass - IRR ButterWorth
% Specifity order = 6
% Fs = 256
% Fc1 = 1.8
% Fc1 = 1.8
sinalFiltrado = filtfilt(SOS1, G1, sinal);

[pxxf, fwf] = pwelch(sinalFiltrado, 1536, fs);

% PSD do sinal filtrado
figure(5);
plot(fwf, pxxf, 'k');
xlim([0 2.5]);
xlabel('Frequência (Hz)');
ylabel('PSD');
title('Densidade espectral de potência do sinal filtrado - BandPass');
grid on;

% Sinal Filtrado 
figure(6);
plot(t, sinalFiltrado)
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal filtrado (BandPass) no canal POz - Trial 3');
grid on;
