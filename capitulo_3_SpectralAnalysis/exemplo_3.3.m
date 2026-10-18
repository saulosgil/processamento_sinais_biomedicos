% Example 3.3 
%
% Plot the magnitude and phase spectrum of a waveform consisting of a single sine wave and white noise 
% with an SNR of −7 dB. Make the waveform 1024-points long. Use MATLAB’s fft routine and plot only the 
% valid points.

% Solution 
%
% Specify N (1024) and fs (1000 Hz) and generate the noisy waveform using sig _ noise. Generate a frequency 
% vector for use in plotting. This vector should be the same length as the signal and ranges to fs 
% (recall that when m = N, f = fs , from Equations 3.5 and 3.6). Take the FT using fft and calculate the 
% magnitude and phase from the complex output. Unwrap the phase and convert it into deg. Plot only the valid 
% points from 1 to N/2. Since the first element produced by the fft routine is the DC component, the frequency
% vector should range between 0 and N–1 if the DC term is included in the spectral plot.

clear all;
close all;
clc;

addpath(genpath("funcoes_auxiliares_and_data_tests"));   % torna sig_noise visível

% Example 3.3 Determine spectrum of a noisy waveform
N = 1024;                                   % Data length
N2 = 511;                                   % Valid spectral points
fs = 1000;                                  % Sample frequency (assumed by sig_noise)
[x,t] = sig_noise(250,-7,N);                % Generate signal (250 Hz sine plus white noise)
X = fft(x);                                 % Calculate FFT
X_mag = abs(X);                             % Compute magnitude spectrum    
Phase = unwrap(angle(X));                   % Phase spectrum unwraped
Phase = Phase*360/(2*pi);                   % Convert phase to deg
f = (0:N-1)*fs/N;                           % Frequency vector

subplot(2,1,1);                 
plot(f(1:N2),X_mag(1:N2),'k');              % Plot magnitude spectrum
title('Magnitude spectrum');
xlabel('Frequency (Hz)');
ylabel('|X(f)|');


subplot(2,1,2);
plot(f(1:N2),Phase(1:N2),'k');              % Plot phase spectrumy
title('Phase spectrum');
xlabel('Frequency (Hz)');
ylabel('Phase(º)');

% Resultados
%
% A presença de uma onda senoidal de 250 Hz é facilmente identificável, mesmo estando “enterrada” 
% em um nível considerável de ruído.

% Os gráficos espectrais de magnitude e fase produzidos pelo programa no Exemplo 3.3. 
% O pico em 250 Hz é evidente. A frequência de amostragem desses dados é de 1 kHz; 
% portanto, o espectro é válido até a frequência de Nyquist, fs/2 (500 Hz). Consequentemente, 
% apenas a primeira metade do vetor produzido pela FFT é representada graficamente (SNR = −7 dB, N = 1024.)

% O sinal analisado
% A função sig_noise(250, -7, N) gera uma senoide de 250 Hz somada a ruído branco, com relação sinal-ruído de −7 dB. Um valor negativo significa que o ruído tem mais potência que o sinal: cerca de 5 vezes mais nesse caso. Se plotar x no tempo, veria só um sinal irregular, e seria praticamente impossível perceber a senoide a olho nu.

% Espectro de magnitude
% Mesmo com o ruído mais forte que o sinal, o pico em 250 Hz se destaca claramente, com altura de ~300, contra um piso de ruído que fica em geral abaixo de 50 a 70. Isso acontece por causa da forma como cada componente distribui sua energia no espectro:
%     - A senoide tem toda a energia concentrada em uma única frequência. Na FFT, ela se acumula em um único ponto e forma um pico alto e estreito.
%     - O ruído branco tem energia espalhada igualmente por todas as frequências. Sua potência total é maior, mas fica dividida entre os 512 pontos do espectro, então cada ponto recebe só uma pequena parte.
% O resultado é que a senoide, fraca no tempo, se torna dominante no ponto exato da sua frequência. Esse “ganho” cresce com o número de amostras: com N maior, o pico fica mais alto em relação ao piso, porque a energia da senoide continua concentrada enquanto a do ruído se espalha por mais pontos.
% O piso de ruído é irregular e aleatório, com picos isolados que podem chegar a 70 ou 90, como em ~45 Hz no livro e ~80 Hz no seu. Eles não representam nenhum componente real: são flutuações estatísticas do ruído, e mudam a cada execução. Na prática, para decidir se um pico é real, compara-se a sua altura com a distribuição do piso, ou se repete a análise em vários trechos e tira-se a média (como no método de Welch que será visto mais à frente do capitulo 3).
% Outros detalhes:
%     - O eixo vai até 500 Hz, que é a frequência de Nyquist (fs/2 = 1000/2). Por isso o código plota só os primeiros N2 = 511 pontos.
%     - A resolução em frequência é fs/N = 1000/1024 ≈ 0,98 Hz, que é o espaçamento entre pontos do espectro.
%     - A escala absoluta (~300) vem da FFT sem normalização, que cresce com N, como no espectro do EEG que você analisou antes.

% Espectro de fase
% Aqui não há um padrão interpretável. Em quase todas as frequências o espectro contém apenas ruído, cuja fase é aleatória em cada ponto. O unwrap acumula esses saltos aleatórios e produz uma curva que vagueia sem direção definida, subindo até ~4000° no livro e descendo até ~−2900° no seu. A diferença entre as duas figuras confirma justamente que essa curva é aleatória.
% A única frequência em que a fase tem significado é 250 Hz, onde ela indica a fase inicial da senoide. Mas esse valor fica perdido no meio da curva acumulada pelo unwrap, e não dá para lê-lo nesse gráfico. Para obtê-lo, seria preciso consultar o valor diretamente, com algo como angle(X(k))*180/pi no índice correspondente a 250 Hz, sem unwrap.

% A mensagem principal
% O exemplo mostra uma das aplicações mais importantes da análise de Fourier: detectar componentes periódicos que estão invisíveis no tempo. Em biossinais isso é muito comum. Um ritmo alfa fraco no EEG, a interferência de 60 Hz da rede elétrica ou respostas evocadas em estado estacionário (SSVEP), em que o cérebro responde na mesma frequência de um estímulo que pisca, podem ficar escondidos sob ruído no tempo e aparecer como picos claros no espectro de magnitude. A figura também reforça que, para sinais ruidosos, a magnitude é a parte informativa do espectro, e a fase só tem significado nas frequências em que existe um componente real.
