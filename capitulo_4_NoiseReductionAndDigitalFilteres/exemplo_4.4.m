%  Example 4.4 Find the actual magnitude spectrum of the rectangular window
% Use the Fourier transform to find the magnitude spectrum of the rectangular window filter 
% given by Equation 4.23 for two different lengths: L = 18 and L = 66. Use a cutoff frequency of 
% 300 Hz assuming a sampling frequency of 1 kHz (i.e., a relative frequency of 0.3).

% Solution
% First, generate the filter’s impulse response, b[k], using Equation 4.19. This can be done by 
% generating the negative half of b, adding b[0] = 2fc, then adding the positive half of b by flipping 
% the negative half. Since the requested L is even, the program rounds up to produce 19 and 67 
% coefficients. After calculating the impulse response, we find the spectrum by taking the Fourier 
% transform of the response and plot the magnitude spectrum.
%
close all; 
clear all;
clc;

N = 256;                                % Numberf of points for plotting
fs = 1000;                              % Sampling freqeuncy
f = (1:N)*fs/N;                         % Frequency vector for plotting
fc = 300/fs;                            % Cutoff frequency (normalized)
L = [18 66];                            % Filter lengths requested

for m = 1:2
    k  = -floor(L(m)/2):-1;             % Construct k for negative b[k];
    b = sin(2*pi*fc*k)./(pi*k);         % Construct negative b[k]
    b =  [b 2*fc, fliplr(b)];           % Rest of b
    H = fft(b,N);                       % Calculate spectrum 

    subplot(1,2,m);                     % Plot magnitude spectrum
    plot(f(1:N/2),abs(H(1:N/2)),'k');
    text(400, 1.0,['Length = ', num2str(length(b))],'FontSize',14);
    xlabel('Frequency (Hz)','FontSize',14);
    axis([0 600 0 1.2]);
    ylabel('B(f)','FontSize',14);
    disp(length(b))
end

% O que as figuras mostram

% As duas curvas tentam reproduzir o retângulo ideal, mas nenhuma consegue com perfeição:

% Banda de transição. A queda de 1 para 0 não é vertical. Com 19 coeficientes, ela vai de ~260 a ~335 Hz. 
% Com 67, fica muito mais íngreme, entre ~290 e ~310 Hz. Quanto mais coeficientes, mais o filtro se aproxima 
% do ideal.

% Ondulações (ripple). Aparecem na banda de passagem (o ganho oscila em torno de 1) e na banda de rejeição 
% (pequenos lóbulos que deixam passar um pouco das frequências acima de 300 Hz). Com 67 coeficientes, as 
% ondulações ficam mais estreitas e concentradas perto do corte.

% O pico perto do corte não diminui. O maior overshoot fica em ~1,09 nas duas figuras, ou seja, 
% cerca de 9% acima de 1, por mais coeficientes que se use. É o fenômeno de Gibbs, o mesmo que você viu na 
% série de Fourier da onda quadrada. A causa também é a mesma: truncar abruptamente. Lá, as frequências da 
% série eram cortadas; aqui, os coeficientes da sinc.