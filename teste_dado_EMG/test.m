close all;
clear all;
clc;

hpf = lerHPF('Sub-C005_tempo-0_tarefa-cvim_braço-nt_medida-torque_Plot_and_Store_Rep_1.0.hpf');
disp(hpf.xmlCanais{1})          % metadados: nomes, unidades, frequência de amostragem
torque = hpf.dados{1};          % escolha o canal certo
fs = hpf.fs(2);                 % confira no XML
t = (0:numel(torque)-1)'/fs;

figure;
plot(t, torque)
xlabel('Tempo (s)');
ylabel('Torque');

% Seleciona apenas o intervalo entre 140 e 160 s
intervalo = t >= 140 & t <= 160;
t_trecho = t(intervalo);
torque_trecho = torque(intervalo);

% Apenas o trecho de 140 a 160 s
figure;
plot(t_trecho, torque_trecho)
xlabel('Tempo (s)');
ylabel('Torque');
title('Trecho de 140 a 160 s');

% Espectro de potência do trecho selecionado
[pxx, fw] = pwelch(torque_trecho, 256, fs);

% PSD do sinal original
figure(2);
plot(fw, pxx, 'k');
xlabel('Frequência (Hz)');
ylabel('PSD');
title('Densidade espectral de potência do sinal original');
grid on;