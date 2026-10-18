% Example 3.4 
% Generate a 0.5-s symmetrical triangle wave assuming a sample frequency of 100 Hz so that N = 50 points. 
% Calculate and plot the magnitude spectrum. Zero pad the signal so that the period is extended to 2 and 8 s 
% and calculate and plot the new magnitude spectra. Limit the spectral plot to a range of 0–10 Hz. 

% Solution 
% First, generate the symmetrical triangle waveform. Since fs = 100 Hz, the signal should be padded by 150 
% and 750 additional samples to extend the 0.5-s signal to 2.0 and 8.0 s. Calculate the appro priate time 
% and frequency vectors based on the padded length. Calculate the spectrum using fft and take the absolute 
% value of the complex output to get the magnitude spectra. Use a loop to calculate and plot the spectra of 
% the three signals. As in Example 3.3, the frequency vector should range from 0 to N–1 since the DC term is 
% included in the spectral plots.

% Example 3.4 Calculate the magnitude spectrum of a simple waveform with % and without zero padding. 
%
clear all;
close all;
clc;

fs = 100;                           % Sample frequencies
N1 = [0 150 750];                   % Padding added to signal
x =[(0:25) (24: - 1:0)];            % Generate basic test signal 
for k = 1:3
    x1 = [x zeros(1,N1(k))];        % Zero pad signal
    N = length(x1);                 % Data length
    t = (1:N)/fs;                   % Time vector for plotting
    f = (0:N-1)*fs/N;               % Frequency vector for plotting

    subplot(3,2,k*2-1);             
    plot(t,x1,'k');                 % Plot test signal
    xlim([0 t(end)]);
    title('Time domain');
    xlabel('Time (s)');
    ylabel('x(t)');

    subplot(3,2,k*2);
    X1 = abs(fft(x1));              % Calculate the magnitude spectrum
    plot(f, X1,'.k');               % Plot magnitude spectrum
    axis([0 10 0 max(X1)*1.2]);
    title('Frequency domain');
    xlabel('Frequency (Hz)');
    ylabel('|X(f)|');
end

% Resultados

% Os gráficos no domínio do tempo e de espectro de magnitude são apresentados na Figura. 
% Todos os gráficos espectrais possuem o mesmo formato, porém, com a técnica de preenchimento com 
% zeros (*zero padding*), os pontos apresentam um espaçamento menor. O simples acréscimo de zeros 
% ao sinal original resulta em uma curva com aspecto visual mais refinado — o que explica a popularidade 
% do *zero padding*, ainda que essa técnica não gere nenhuma informação adicional. Neste exemplo, o sinal 
% foi preenchido com zeros antes da chamada da rotina FFT, uma vez que se desejava obter os gráficos no 
% domínio do tempo dos sinais assim processados. Alternativamente, o segundo argumento da rotina FFT poderia 
% ter sido utilizado para realizar esse preenchimento, conforme descrito anteriormente. Além disso, o espectro
% completo foi plotado, incluindo os pontos redundantes; contudo, tais pontos foram removidos por meio da 
% limitação do eixo x, utilizando a função `xlim`.

% A figura mostra o efeito do zero padding: acrescentar zeros ao final de um sinal antes de calcular a FFT. 
% O mesmo pulso triangular, com duração de 0,5 s e pico de 25, é analisado três vezes, com comprimentos 
% totais diferentes.

% O que acontece no espectro
% O espaçamento entre os pontos da FFT é sempre 1/T. Ao acrescentar zeros, T aumenta e os pontos ficam mais 
% próximos uns dos outros.
% No painel superior, sem padding, só há pontos a cada 2 Hz: 0, 2, 4, 6, 8 e 10 Hz. Com tão poucos pontos, 
% é difícil perceber o formato do espectro. Na prática, a FFT está tratando o triângulo como se ele se 
% repetisse a cada 0,5 s, e esses pontos são os harmônicos dessa onda periódica.
% No painel do meio, com pontos a cada 0,5 Hz, o formato já começa a aparecer: uma queda suave até ~4 Hz e 
% um pequeno lóbulo entre 4 e 8 Hz.
% No painel inferior, com pontos a cada 0,125 Hz, a curva fica praticamente contínua e o formato do espectro 
% aparece por completo.
% O mais importante é que os três espectros descrevem a mesma curva. Os pontos do painel superior caem 
% exatamente sobre a curva do painel inferior. O padding não altera o espectro: só amostra a mesma curva 
% com mais densidade, ou seja, faz uma interpolação. É isso que o texto do livro diz: o padding produz uma 
% curva com aparência melhor, mas não gera informação nova.

% O formato do espectro
% A curva tem formato de sinc². Isso acontece porque um triângulo é a convolução de dois pulsos retangulares, 
% e convolução no tempo equivale a multiplicação na frequência: sinc × sinc = sinc².
% Os zeros caem em 4 e 8 Hz. Cada pulso retangular que forma o triângulo dura 0,25 s, e os zeros da sinc 
% ficam em múltiplos de 1/0,25 = 4 Hz, a mesma lógica do pulso retangular que você analisou antes. 
% Comparado com o espectro do pulso retangular, os lóbulos laterais do triângulo são bem menores e 
% decaem mais rápido (com 1/f² em vez de 1/f). Isso acontece porque o triângulo não tem descontinuidades: 
% ele tem apenas “quinas”, e transições mais suaves no tempo exigem menos energia em frequências altas.

% Por que a magnitude é a mesma nos três casos
% O valor em 0 Hz é ~625 nos três painéis. Em 0 Hz, a FFT é simplesmente a soma de todas as amostras. 
% Os zeros acrescentados não contribuem para essa soma, então ela não muda: 625 ≈ área do 
% triângulo × fs = (0,5 × 25 / 2) × 100.

% O que o zero padding faz e o que não faz
% O padding melhora a aparência e a leitura do espectro. Fica mais fácil localizar picos, zeros e lóbulos, 
% e estimar com mais precisão a frequência de um pico que caia entre dois pontos.
% O padding não melhora a resolução real, ou seja, a capacidade de separar duas frequências próximas. 
% Essa capacidade depende da duração do sinal medido, e não do número de pontos da FFT. Se duas oscilações 
% estiverem a 0,5 Hz de distância e você registrou apenas 0,5 s de sinal, os picos vão aparecer fundidos, 
% por mais zeros que você acrescente. A única forma de separá-los é registrar um trecho mais longo.

% # IMPORTANTE #
% Como o livro menciona, não é necessário acrescentar os zeros manualmente: fft(x, 800) faz o padding 
% automaticamente até 800 pontos. O livro fez manualmente só para poder plotar os sinais com padding no 
% domínio do tempo.
