% Exemplo 2.10
%
% Compare the EEG signal shown in Figure 2.15 with sinusoids ranging in frequencies between 
% 1.0 and 25 Hz. The sinusoidal frequencies should be in 1.0 Hz increments.

% Solution
%
% Load the ECG signal found as vector eeg in file eeg_data.mat. Use a loop to generate a 
% series of sine waves from 0.25 to 25 Hz. Since we do not know the best frequencies to use in the
% correlation, we will increment the sine wave frequency in small intervals of 0.25 Hz. (Cosine 
% waves would work just as well since the cross-correlation covers all possible phase shifts.) 
% Cross-correlate these sine waves with the EEG signal and find the maximum cross-correlation. 
% Plot this maximum correlation as a function of the sine wave frequency.

% Example 2.10 Camparison of an EEG signal with sinusoids
%
clear all;
close all;
clc;

caminho = 'funcoes_auxiliares_and_data_tests\Chapter 2\eeg_data.mat'
load(caminho);

fs = 50;                            % frequência de amostragem (Hz)
t = (1:length(eeg))/fs;             % vetor de tempo

for i = 1:100
    f(i) = 0.25 * i;                % 0,25 a 6,25 Hz
    x = sin(2 * pi * f(i) * t);     % senoide de teste
    rxy = axcor(eeg, x);            % correlação cruzada
    rmax(i) = max(rxy);             % maior correlação entre todos os deslocamentos
end

plot(f, rmax, 'k'); 
xlabel('Frequency(Hz)');
xticks(0:2.5:25);
ylabel('Maximum correlation');
grid on;

% Resultado

% O resultado das múltiplas correlações cruzadas é mostrado na Figura 2.16, e 
% surge uma estrutura interessante. Algumas frequências apresentam correlação 
% muito maior entre a senoide e o EEG. Um pico particularmente forte é observado 
% na região de 2.5–7.5 Hz, indicando a presença de um padrão oscilatório conhecido 
% como onda alfa. 
% 
% A transformada de Fourier é um método mais eficiente para obter 
% a mesma informação, como mostrado no Capítulo 3.

% O que o gráfico mostra

% A atividade está concentrada abaixo de ~10 Hz. Entre 0 e 9,5 Hz a correlação fica quase toda ]
% acima de 0,13, com picos entre 0,16 e 0,24. Logo depois de 9,5 Hz há uma queda abrupta, 
% e de 10 a 25 Hz a curva oscila em um patamar baixo, entre 0,05 e 0,11, sem nenhum pico relevante. 
% Em outras palavras, esse EEG é dominado por ritmos lentos, e há pouca atividade beta (13–30 Hz).

% Dentro da faixa baixa aparecem três regiões de destaque:
    % 2,5 a 3,75 Hz (delta/teta baixo): dois picos de ~0,23. Atividade lenta pode refletir sonolência, 
    %   mas também artefatos como movimentos oculares e piscadas, que costumam contaminar essa faixa.
    % 6,25 a 6,75 Hz (fronteira teta/alfa): os picos mais altos do gráfico, ~0,24.
    % 8 a 9 Hz (alfa): picos secundários de ~0,165, que correspondem à onda alfa citada no livro.