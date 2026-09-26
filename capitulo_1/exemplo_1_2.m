% Exemplo 1.2  Gera uma onda senoidal discreta de 2 Hz usando MATLAB.
% Assume um tempo de amostragem de 0.01 seg. e use pontos suficientes para tornar
% a onda senoidal 1 seg. de comprimento; i.e., o tempo total, TT deve ser 1 seg.

clear all; % limpa todas variaveis do workspace
close all; % fecha todas as figuras abertas

Ts = .01;               % Define o tempo de amostragem    
TT = 1;                 % Define o tempo total
f = 2;                  % Define a frequência
t = 0:Ts:TT;            % Gera o vetor de tempo
x = sin(2*pi*f*t);      % Gera o sinal senoidal desejado

plot(t,x,'.k');         % Plota o sinal
xlabel('Time (sec)');   % Rotula o eixo x
ylabel('x(t)');         % Rotula o eixo y