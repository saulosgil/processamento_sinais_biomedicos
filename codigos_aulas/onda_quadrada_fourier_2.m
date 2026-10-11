% Estima coeficientes Transformada de Fourier (TF) e Série de Fourier (SR)
% 12/10/2026
% Saulo Gil

clear all;
close all;
clc;

T = 2;                      % periodo (sec)
fs  = 500;                  % frequecnia de amostragem
f1 = 1/T;                   % frequeencia fudamental
N = round(T * fs);          % numero de pontos
t = (1:N) / fs;             % vetor de tempo

% Gerando o pulso retangular no tempo
x = zeros(1, N);
PW = round(0.05 * fs);
x(1:PW) = 1;

figure(1);
plot(t, x, 'linewidth', 2);
xlabel('Tempo (sec)');
ylabel('Amplitude (a.u.)');
title('Sinal no Tempo');
grid on;

% Coeficientes da SF e definição da TF

for m = 1:round(N/2)
    f(m) = m * f1;
    C(m) = mean(x.*exp(-1i * 2 * pi * m * f1 * t));
    Cmag(m) = sqrt(real(C(m))^2 + imag(C(m))^2);
    Xfreq(m) = T*Cmag(m);
end

figure(2);
plot(f, Cmag, '.r');
xlabel('Frequencia (Hz)');
ylabel('Cmag (a.u.)');
title('Espectro de Magnitude');
grid on;

figure(3);
plot(f, Xfreq, '.g');
xlabel('Frequencia (Hz)');
ylabel('Xfreq (f)');
title('Transformada de Fourier');
grid on;


% As três figuras mostram um pulso retangular no tempo e o seu conteúdo em frequência, 
% calculado de duas formas: pelos coeficientes da série de Fourier e pela transformada de Fourier.

% Figura 1: o sinal no tempo

% É um pulso de amplitude 1 e duração de 0,05 s (50 ms), seguido de zeros até completar 2 s. O registro total de 2 s define o período considerado na análise, T = 2 s. Pelo eixo de frequência, que vai até 250 Hz, a frequência de amostragem é de 500 Hz.

% Figuras 2 e 3: o espectro do pulso

% As duas curvas têm o mesmo formato, o da função sinc em módulo, |sen(x)/x|:

% um lóbulo principal grande, centrado em 0 Hz;
% uma sequência de lóbulos laterais cada vez menores;
% zeros em 20, 40, 60, 80 Hz e assim por diante.

% Por que os zeros caem em múltiplos de 20 Hz? A posição dos zeros é o inverso da duração do pulso: 1/0,05 s = 20 Hz. 
% Uma senoide de 20 Hz completa exatamente um ciclo dentro do pulso, e por isso sua integral sobre ele dá zero. O mesmo 
% vale para 40 Hz (dois ciclos), 60 Hz (três ciclos) etc.

% A relação inversa entre tempo e frequência. Essa é a ideia central da figura. Quanto mais curto o pulso, mais largo o 
% lóbulo principal e mais espalhada a energia pelas frequências. Um pulso de 10 ms teria o primeiro zero em 100 Hz. Um de 200 ms,
% em 5 Hz. No limite, um impulso instantâneo tem espectro completamente plano: contém todas as frequências com a mesma amplitude. 
% Por isso eventos rápidos (picos, bordas, artefatos) espalham energia por uma faixa larga do espectro.

% Os lóbulos laterais. O primeiro lóbulo lateral (~30 Hz) tem cerca de 22% da altura do principal (0,0055 / 0,025), e os 
% % seguintes decaem aproximadamente com 1/f. Esse decaimento lento é consequência das bordas abruptas do pulso. É o mesmo 
% % fenômeno da série de Fourier da onda quadrada e da rampa: descontinuidades exigem muitas frequências altas.

% O módulo esconde o sinal da função. A sinc real alterna entre positivo e negativo de um lóbulo para o outro. Como o gráfico 
% mostra o valor absoluto, todos os lóbulos aparecem positivos. Essa alternância está guardada no espectro de fase, como saltos 
% de 180° em cada zero.

% Por que a curva aparece pontilhada? Cada ponto é uma frequência calculada. O espaçamento entre eles é a resolução em frequência, 
% 1/T = 1/2 s = 0,5 Hz. Um registro mais longo, com mais zeros depois do pulso, deixaria os pontos mais próximos, mas não mudaria o 
% formato da curva.

% A diferença entre as Figuras 2 e 3

% O formato é idêntico. Só a escala muda, e a Figura 3 tem exatamente o dobro dos valores da Figura 2 (0,05 contra 0,025 em 0 Hz). 
% A diferença está na normalização:

% Figura 2 (série de Fourier): os coeficientes são divididos pelo período. O valor em 0 Hz é a média do sinal ao longo dos 2 s: 
% área / T = 0,05 / 2 = 0,025.

% Figura 3 (transformada de Fourier): não divide pelo período. O valor em 0 Hz é a área do pulso: 1 × 0,05 s = 0,05.

% As duas grandezas se relacionam por X(f) = T · cₙ. Com T = 2 s, a transformada é o dobro dos coeficientes. Em termos de 
% interpretação, os coeficientes da série descrevem a amplitude de cada harmônico de um sinal periódico, enquanto a transformada 
% descreve a densidade espectral de um sinal não periódico, com unidade de amplitude × tempo.

% A mensagem principal

% O par pulso retangular ↔ sinc ilustra o princípio mais importante da análise de Fourier: o que é estreito em um domínio é largo 
% no outro. Esse princípio tem consequência prática direta no seu trabalho com EEG. Uma janela de análise curta tem boa resolução 
% no tempo, mas espectro “borrado”. Uma janela longa tem resolução fina em frequência, mas perde a localização temporal dos eventos. 
% Os lóbulos laterais da sinc também explicam o vazamento espectral (leakage) visto antes: cortar um sinal com uma janela retangular 
% equivale a multiplicá-lo por esse pulso, e os lóbulos laterais espalham energia para frequências vizinhas. Janelas suaves, como 
% Hanning e Hamming, existem justamente para reduzir esses lóbulos.