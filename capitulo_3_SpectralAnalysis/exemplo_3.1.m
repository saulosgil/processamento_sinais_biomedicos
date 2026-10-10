% Example 3.1
%
% Use Fourier series analysis to generate the magnitude and phase plot of the ECG signal origi
% nally shown in Figure 2.16 and repeated below. For these data, fs = 50 Hz.

% Solução

% A análise da série de Fourier requer que a forma de onda seja periódica; portanto, 
% devemos assumir que o sinal de EEG é periódico. Em outras palavras, assume-se que o 
% segmento de EEG corresponde a um ciclo de um sinal periódico. Obviamente, isso é 
% extremamente improvável, mas não importa, pois não sabemos nada sobre o sinal antes ou 
% depois do segmento no computador; assim, qualquer suposição sobre esse sinal 
% desconhecido é aceitável. 

% Feita essa suposição, a implementação das Equações 3.3 e 3.4 no MATLAB torna-se simples. 
% Em seguida, utilizamos as Equações 3.11 e 3.12 para converter os componentes seno e cosseno 
% em gráficos de magnitude e fase. 

% Definimos M = N/2, pois esse valor corresponde à frequência máxima válida encontrada no 
% sinal digitalizado (Equação 3.7). (Como essa frequência é fs/2 e fs = 50 Hz, prevemos uma 
% frequência máxima de 25 Hz.)

% Example 3.1 Use Fourier series analysis to generate the magnitude and
% phase plots of the ECG signal.
%
clear all;
close all;
clc;

fs = 50; 

caminho = 'funcoes_auxiliares_and_data_tests\Chapter 2\eeg_data.mat';
load(caminho); 

N = length(eeg);
Tt = N/fs; 
f1 = 1/Tt; 
t = (1:N)/fs; 

for m = 1:round(N/2)
    f(m) = m * f1;
    a = sum(eeg.*cos(2 * pi * f(m) * t));
    b = sum(eeg.*sin(2 * pi * f(m) * t));
    X_mag(m) = sqrt(a^2 + b^2);
    X_phase(m) = -atan2(b,a);
end

X_phase = unwrap(X_phase); 
X_phase = X_phase * 360 / (2 * pi); 

subplot(2,1,1);
plot(f,X_mag,'k'); 
ylabel('|X(m)|');
xlabel('Frequency (Hz)');

subplot(2,1,2);
plot(f,X_phase,'k'); 
ylabel('Phase(º)');
xlabel('Frequency (Hz)');
    
% Espectro de magnitude (painel superior)

% Indica quanto de cada frequência existe no sinal.
%     -   6 a 9 Hz: a região dominante, com o pico mais alto em ~6,3 Hz e picos fortes até ~8,5 Hz. Corresponde à transição teta/alfa e à banda alfa.
%     -   1,5 a 2,5 Hz: um grupo de picos menores na faixa delta.
%     -   16 a 20 Hz: um grupo bem definido na banda beta. Esse grupo era quase invisível no método de correlação e aparece com clareza na FFT, que é mais sensível e precisa.
%     -   Acima de 21 Hz: pouca atividade.

%     Por que a curva é tão serrilhada? Cada ponto da FFT de um único trecho de EEG é uma 
%     estimativa com muita variância. O sinal é aleatório e não estacionário, e nenhum 
%     janelamento ou média foi aplicado. Por isso se olha mais para o formato geral das 
%     regiões do que para picos individuais. Métodos como o de Welch (pwelch) dividem o 
%     sinal em trechos, calculam o espectro de cada um e tiram a média, produzindo uma curva 
%     bem mais suave.
% ----------------------------------------------------------------------------------------------
% Espectro de fase (painel inferior)

% Mostra a defasagem de cada componente de frequência. Os valores chegam a −14.000° porque a 
% fase foi desembrulhada (unwrap): em vez de ficar presa entre −180° e +180°, ela vai acumulando 
% voltas completas de 360°.

% Para um sinal como o EEG, a fase sozinha tem pouco valor interpretativo. O EEG se comporta 
% como um processo aleatório, e a fase de cada componente depende de onde o trecho registrado 
% começa. A tendência descendente e aproximadamente linear é, em boa parte, efeito do próprio 
% desembrulhamento sobre fases quase aleatórias, e não um fenômeno fisiológico. 

% Por isso, em análise de EEG, quase sempre se trabalha com a magnitude ou a potência. 
% A fase volta a ser útil quando se comparam sinais entre si, por exemplo na coerência ou 
% na sincronia de fase entre canais, ou em potenciais evocados, onde a resposta está alinhada 
% ao estímulo.