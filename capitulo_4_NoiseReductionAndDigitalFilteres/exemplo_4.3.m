% Example 4.3 Spectra of moving average filters.
% Find the magnitude and phase spectra of a 3-point moving average filter and a 10-point moving 
% average filter, both with equal weights. Plot the two spectra superimposed. Assume a sampling 
% frequency of 250 Hz.

% Solution
% The frequency spectrum of any FIR filter can be determined by taking the Fourier transform 
% of its filter coefficients, which is also the impulse response (Equation 4.15). However, since the 
% impulse responses of many FIR filters can be quite short (only three points for a three-point 
% moving average filter), it is appropriate to pad the response before taking the Fourier transform. 
% In this example, we first construct the two filter impulse responses, then take the Fourier trans
% form using fft, but pad the data out to 256 points to get smooth plots. We then plot the valid 
% points of the magnitude and phase from the complex Fourier transform using a frequency vec
% tor based on the 250 Hz sampling rate.
clear all; 
close all; 
clc;

fs = 250;               % Sampling frequency      

%  Construct moving average filters
b3 = [1 1 1]/3;
b10 = [1 1 1 1 1 1 1 1 1 1]/10;      % Construct filter impulse response
%
B3 = fft(b3,256);                        % Find filter system spectrum
B3_mag = abs(B3);                        % Compute magnitude spectrum
Phase3 = unwrap(angle(B3));              % Compute phase angle, unwrap
Phase3 = Phase3*360/(2*pi);              % Phase angle in degrees
%
B10 = fft(b10,256);                      % Find filter system spectrum
B10_mag = abs(B10);                      % Compute magnitude spectrum
Phase10 = unwrap(angle(B10));            % Compute phase angle, unwrap
Phase10 = Phase10*360/(2*pi);            % Phase angle in degrees

% Plot filter spectrum magnitude and phase
N1 = round(length(B3_mag)/2);            % Number of frequency points to plot
f = (1:N1)*fs/(2*N1);                    % Frequency vector for plotting spectrum

subplot(2,1,1); hold on;
plot(f,B10_mag(1:N1),'k')                % Plot magnitude 
plot(f,B3_mag(1:N1),':k','LineWidth',1)  % Plot magnitude 
xlabel('Frequency (Hz)','FontSize',14); 
ylabel('|B|','FontSize',14);
xlim([0 f(N1)]);

subplot(2,1,2); hold on;
plot(f,Phase10(1:N1),'k');               % Plot phase
plot(f,Phase3(1:N1),':k','LineWidth',1); 
xlabel('Frequency (Hz)','FontSize',14); 
ylabel('Phase (deg)','FontSize',14);
xlim([0 f(N1)]);

% Magnitude (painel superior)

% Os dois filtros são passa-baixas: deixam passar 0 Hz com ganho 1 e atenuam as frequências mais altas. 
% Isso faz sentido, porque fazer a média suaviza o sinal e remove as variações rápidas.

% O filtro de 10 pontos corta bem mais cedo. A magnitude cai pela metade por volta de 11 Hz, enquanto no de 
% 3 pontos isso só acontece por volta de 50 Hz. Quanto mais pontos na média, mais forte a suavização.

% Os zeros caem em múltiplos de fs/L:

%     - 3 pontos: 250/3 ≈ 83 Hz;
%     - 10 pontos: 250/10 = 25, 50, 75 e 100 Hz.

% Nessas frequências, a janela de média contém exatamente um número inteiro de ciclos da senoide, e a média 
% dá zero. Isso é útil na prática: uma média móvel com o comprimento certo elimina uma frequência específica. 
% Por exemplo, uma média de 5 pontos com fs = 300 Hz zeraria os 60 Hz da rede elétrica.

% O formato é de uma sinc. A média móvel é um pulso retangular, e você já viu que o espectro de um retângulo é 
% uma sinc. Por isso aparecem os lóbulos laterais (o primeiro chega a ~0,22 no filtro de 10 pontos). 
% A média móvel é um passa-baixas simples, mas de baixa qualidade: ela não corta de forma limpa, e parte das 
% frequências altas passa pelos lóbulos.

% Fase (painel inferior)

% A fase cai em linha reta. Isso é a fase linear, e corresponde a um atraso constante para todas as frequências:

%     - 3 pontos: atraso de (3−1)/2 = 1 amostra (4 ms);
%     - 10 pontos: atraso de (10−1)/2 = 4,5 amostras (18 ms).

% A reta do filtro de 10 pontos é mais inclinada porque o atraso é maior. Fase linear é uma qualidade 
% importante: todas as frequências chegam atrasadas pelo mesmo tempo, então o formato do sinal não é distorcido, 
% apenas deslocado. Filtros FIR com coeficientes simétricos sempre têm essa propriedade, ao contrário de 
% filtros IIR como o do exemplo anterior.

% Os saltos verticais de +180° aparecem exatamente nos zeros da magnitude (25, 50, 75, 100 Hz e 83 Hz). 
% Neles, a sinc muda de sinal, de positivo para negativo, e trocar o sinal equivale a somar 180° na fase. 
% O unwrap só corrige saltos de 360°, por isso esses de 180° continuam no gráfico. Eles não indicam nenhum 
% problema: entre os saltos, a fase continua linear.

% Em resumo

% Aumentar o tamanho da média móvel tem um compromisso claro:

%     - mais pontos: mais suavização, corte em frequência mais baixa e mais atraso;
%     - menos pontos: preserva mais detalhes do sinal, mas remove menos ruído.