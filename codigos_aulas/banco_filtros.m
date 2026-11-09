% Atividade filtragem digital de sinais
clear all;
close all;
clc;

addpath(genpath("atividades_praticas"));   % torna broadband1.mat visível

load('S03.mat'); % Carrega os filtros projetados
load('CoefButterTest.mat')

[Nsamples, Nchannels,trials] = size(stimuliData); % Obtém o número de amostras e canais do sinal

t = 0:(1/fs):(Nsamples-1)*(1/fs);  % Vetor de tempo
sinal10 = stimuliData(:, 5, 3) % POz, Trial 3, estimulo 10Hz

% Visualização do canal POz para o trial 3
figure;
plot(t, sinal10);
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal no canal POz - Trial 3');
grid on;

% fdatool % Abre a toolbox de projeto de filtros OBS. Rodar no terminal do
% MATLAB para projetar os filtros. Depois, salvar os filtros
% (SOS e G; CoefButterTest.mat)

% Banco de filtros
sinalFiL10 = filtfilt(SOS10, G10, sinal10);
sinalFiL11 = filtfilt(SOS11, G11, sinal10);
sinalFiL12 = filtfilt(SOS12, G12, sinal10);
sinalFiL13 = filtfilt(SOS13, G13, sinal10);

figure;
subplot(2,2,1)
plot(t, sinalFiL10, 'DisplayName', '10 Hz');
xlabel('Tempo (s)'); ylabel('Amplitude');
title('Sinais filtrados - Trial 3'); grid on; legend;

subplot(2,2,2)
plot(t, sinalFiL11, 'DisplayName', '11 Hz');
xlabel('Tempo (s)'); ylabel('Amplitude');
title('Sinais filtrados - Trial 3'); grid on; legend;

subplot(2,2,3)
plot(t, sinalFiL12, 'DisplayName', '12 Hz');
xlabel('Tempo (s)'); ylabel('Amplitude');
title('Sinais filtrados - Trial 3'); grid on; legend;

subplot(2,2,4)
plot(t, sinalFiL13, 'DisplayName', '13 Hz');
xlabel('Tempo (s)'); ylabel('Amplitude');
title('Sinais filtrados - Trial 3'); grid on; legend;

% Comparação espectral dos sinais filtrados
figure;
periodogram([sinalFiL10, sinalFiL11, sinalFiL12, sinalFiL13], [], [], fs);
legend('10 Hz', '11 Hz', '12 Hz', '13 Hz');
title('Espectro dos sinais filtrados');

% Potência de cada sinal filtrado - P10
potencia10 = mean(sinalFiL10.^2);
potencia11 = mean(sinalFiL11.^2);
potencia12 = mean(sinalFiL12.^2);
potencia13 = mean(sinalFiL13.^2);

potencias = [potencia10, potencia11, potencia12, potencia13];
potencias