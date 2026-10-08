% Esse é o Exemplo 1.4 do livro do Semmlow (Biosignal and Medical Image Processing). 
% Ele compara o ruído de quantização medido na prática com o valor teórico previsto 
% pela Eq. 1.8, ou seja, verifica se a fórmula q²/12 descreve bem o erro introduzido 
% ao digitalizar um sinal com poucos bits.

% Example 1.4 Avalia a equação de quantização utilizando dados simulados.

f = 4;    % frequencia
N = 1000;   % Número de pontos
Ts = 0.002;   % Ts
bits = 6;   % Nível de quantização
t = (0:N-1)*Ts;   % Vetor usado para gerar onde senoidal(1-cycle)
signal_in = sin(2*pi*f*t); % Gera o sinal

signal_out = quantization(signal_in, bits); % Quantiza o sinal
noise_signal = signal_out - signal_in; % Determina o erro de quantização
q_noise = var(noise_signal); % Variancia do erro de quantização
q = 1/(2^bits - 1);  % Calcula oNível de quantização (Eq. 1.6)
theoretical = (q^2)/12; 

disp(' Quantization Noise')
disp('Bits Emperical Theoretical') % quantization error (Eq. 1.8)
out = sprintf('%2d %5e %5e', bits, q_noise, theoretical); % Format output
disp(out)