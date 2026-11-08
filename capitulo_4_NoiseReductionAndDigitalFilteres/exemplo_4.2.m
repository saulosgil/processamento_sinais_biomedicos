% Example 4.2
% A system is defined by the digital transfer function below. Plot its frequency spectrum: mag
% nitude in dB and phase in degrees. Also plot the time-domain response of the system to an 
% impulse input (i.e., the system’s impulse response). Assume a sampling frequency of fs = 1 kHz 
% and make N = 512. These are both optional: fs is provided so that the frequency curve can be 
% plotted in Hz instead of relative frequency, and N is chosen to be large enough to produce a 
% smooth frequency curve. (In this example, it is actually much larger than needed for a smooth 
% curve: N = 128 would be sufficient).

% Solution
% 4.4 Finite Impulse Response Filters
% In filter design, we usually work from a desired frequency spectrum, determining the a and b 
% weights that generate a digital filter to approximate the desired spectrum. In this example, we 
% work the other way around: we are given the coefficients in terms of the digital transfer func
% tion and are asked to find the spectrum. From the transfer function, we note that a[0] = 1 (as 
% required), a[1] = 0.2, a[2] = 0.8, b[0] = 0.2, and b[1] = 0.5. We apply Equation 4.13 to find H[m], 
% then convert to frequency in Hz by multiplying m by fs
% /N. To find the impulse response, we 
% generate an impulse function (1.0 followed by zeros) and determine the output by implementing 
% Equation 4.12 through the MATLAB filter routine.
close all; 
clear all;
clc;

fs = 1000;					            % Sampling frequency
N = 512;					            % Number of points

%  Define a and b coefficients based on H(z) - Vide Equation
a = [1 -.2 .8];
b = [.2 .5];

% Plot the Frequency characteristic of H(z) using the fft
H = fft(b,N)./fft(a,N);		            % Compute H(f)
Hm = 20*log10(abs(H));                  % Get magitude in dB
Theta = (angle(H))*360/(2*pi);          % and phase in deg. 
f = (1:N)*fs/N;				            % Frequency vector for plotting

subplot(2,1,1); 
plot(f(1:N/2),Hm(1:N/2),'k');			% Plot and label mag  H(f)
xlabel ('Frequency (Hz)','FontSize',14);
ylabel('|H(z)| (dB)','FontSize',14);
grid on;                                % Plot using grid lines

subplot(2,1,2); 
plot(f(1:N/2),Theta(1:N/2),'k');		% Repeat for phase 
xlabel ('Frequency (Hz)','FontSize',14);
ylabel('Phase (deg)','FontSize',14);
grid on;                                % Plot using grid lines

% Compute the Impulse Response
x = [1, zeros(1,N-1)];		            % Generate an impulse function
y = filter(b,a,x);				        % Apply b and a to impulse using Eq. 4-6 

figure;						            % New figure
t = (1:N)/fs;
plot(t(1:60),y(1:60),'k'); hold on;		% Plot only the first 60 points for clarity
plot([t(1) t(60)],[0 0],'k');           % Plot zero line
xlabel('Time (sec)','FontSize',14);		% Labels
ylabel ('Impulse Response','FontSize',14);

% Figura 1: magnitude e fase

% Magnitude: o filtro tem um pico de ressonância em ~233 Hz, com ganho de +9 dB. Ele amplifica as frequências 
% perto desse valor e atenua as outras: −7 dB nas baixas e −16 dB perto de 500 Hz. Funciona como um passa-banda, 
% que destaca uma faixa estreita de frequências.

% Fase: é praticamente 0° até ~200 Hz e cai rapidamente em torno do pico, chegando a ~−180°. Essa queda brusca 
% de fase justo na frequência de ressonância é típica de sistemas ressonantes.

% O pico vem do denominador. As raízes de z² − 0,2z + 0,8 = 0 (os polos) ficam a um ângulo que corresponde a 
% 232 Hz e bem perto do círculo unitário (raio 0,89). Quanto mais perto do círculo, mais alto e mais estreito 
% o pico.

% Figura 2: resposta ao impulso

% Depois do impulso, o filtro oscila e a oscilação vai diminuindo até sumir, em ~50 ms. É como um sino que, 
% depois de tocado, vibra e se cala aos poucos.

% As duas figuras mostram a mesma coisa vista de lados diferentes:

% A frequência da oscilação é a do pico. Os picos da curva se repetem a cada ~4,3 ms, e 1/0,0043 ≈ 232 Hz.
% A velocidade com que ela some depende do raio dos polos (0,89). Polos mais perto do círculo dariam uma 
% oscilação mais longa e um pico mais estreito na magnitude.
% A oscilação nunca chega exatamente a zero, por causa da realimentação. Daí o nome “resposta ao impulso 
% infinita”.