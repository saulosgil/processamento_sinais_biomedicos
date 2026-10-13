% SF (Série de Foutier) e DFT (Discrete Fourier Transform) para onda triangular
% 13/10/2026
% Saulo Gil

close all;
clear all;
clc;

T = 1;                              % sex
fs = 256;                           % freq de amostragem
N = fs*T;                           % Número de pontos
t = (1:N)*T/N;                      
f = (1:N)*fs/N;

x = [t(1:N/2) zeros(1,N/2)];

figure(1);
plot(t,x);
xlabel('tempo (s)');
ylabel('Ampliture');
grid on;

% Fourier Transform - via fft
Xf = fft(x);
Mag = abs(Xf(2:end))/(N/2);
Phase = angle(Xf(2:end))*360/(2*pi);

figure(2);
plot(f(1:20), Mag(1:20), 'xb');
xlabel('Frequencia (Hz)');
ylabel('|X(f)|');
grid on;

% Calculando os coeficiêntes para série de Fourier
for m = 1:20
    a(m) = (2/N)*sum(x.*(cos(2*pi*m*t)));
    b(m) = (2/N)*sum(x.*(sin(2*pi*m*t)));
    C(m) = sqrt(a(m).^2 + b(m).^2);
    theta(m) = 360/(2*pi)*atan2(-b(m), a(m));
end

disp([a(1:4)' b(1:4)' C(1:4)' theta(1:4)' Mag(1:4)' Phase(1:4)'])

figure(2); hold on;
plot(f(1:20), C(1:20), 'sr');
