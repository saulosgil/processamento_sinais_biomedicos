% Exemplo 2.6
% Generate a 500 point, 2 Hz sine wave and a 4 Hz sine wave of the same length. Make TT = 1 s. 
% Are these two waveforms orthogonal?

% Solution
% Generate the waveforms. 
% Since N = 500 and TT = 1 s, Ts = TT/N = 0.002 s. 
% Apply Equation 2.29 to these waveforms. 
% If the result is near zero, the waveforms are orthogonal.
clear all;
close all;
clc;

% Example 2.6 Evaluate 2 waveform for Orthogonality.
Ts = 0.002;           % intervalo de amostragem 
N = 500;              % Numero de amostras
t = (0:N-1)*Ts;       % vetor de tempo
f1 = 2;               % frequência da primeira onda
f2 = 4;               % frequência da segunda onda
x = sin(2*pi*f1*t);   % primeira onda
y = sin(2*pi*f2*t);   % segunda onda
Corr = sum(x.*y);     % correlação entre as duas ondas
disp(Corr)            % exibe o resultado da correlação (notação científica)
fprintf('Correlação: %.4f\n', Corr)   % exibe em decimal fixo

% Result
% The correlation value produced by application of Equation 2.29 to the two sine waves is -5.7200e-15. 
% This is very close to zero and shows the waveforms are orthogonal. 
% This is expected, since harmonically related sines (or cosines) are known to be orthogonal. 
% This is one of the reasons for % the utility of sinusoids in the Fourier transform as detailed 
% in Chapter 3.