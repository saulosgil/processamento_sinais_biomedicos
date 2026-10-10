% Exemplo 3.2
% Perform a discrete Fourier series analysis on the triangular waveform defined by the equation:

%         t 0   <= t < 0.5 s 
% x(t) = 
%         t 0.5 <= t < 1.0 s

% Reconstruct the waveform using the first five and then the first 10 components. Plot the 
% reconstructed waveforms superimposed on x(t). Make fs = 500 Hz.

% Solução

% Após construir a forma de onda no MATLAB, decomponha-a em 10 componentes usando a abordagem 
% apresentada no Exemplo 3.1. 
% Observe que o tempo total da forma de onda é de 1,0 s. Em seguida, reconstrua as duas formas
% de onda usando a versão digital da Equação 3.15. No entanto, se usarmos a Equação 3.15 para a
% reconstrução, precisaremos normalizar de maneira equivalente à usada nas Equações 3.1 e 3.2, 
% ou seja, por 2/N. Observe que a forma de onda possui um termo CC que deve ser adicionado às 
% reconstruções.

clear all;
close all;
clc;

% Example 3.2 Fourier series decomposition and reconstruction.
%
fs = 500;                               % Sampling frequency
Tt = 1;                                 % Total time
N = Tt*fs;                              % Determine N
f1 = 1/Tt;                              % Fundamental frequency
t = (1:N)/fs;                           % Time vector
x = zeros(1,N);                         % Construct waveform
x(1:N/2) = t(1:N/2);                    % Rampa de 0 a 0,5 s

% Decomposição de Fourier
n_max = 100;                             % Maior harmônico calculado
f = zeros(1,n_max);
X_mag = zeros(1,n_max);
X_phase = zeros(1,n_max);

a0 = 2*mean(x);                         % Calculate a(0)

for m = 1:n_max
   f(m) = m*f1;                         % Sinusoidal frequencies
   a = (2/N)*sum(x.*cos(2*pi*f(m)*t));  % Cosine coeff., Eq. 3.3 
   b = (2/N)*sum(x.*sin(2*pi*f(m)*t));  % Sine coeff., Eq. 3.4
   X_mag(m) = sqrt(a^2 + b^2);          % Magnitude spectrum Eq. 3.6
   X_phase(m) = -atan2(b,a);            % Phase spectrum, Eq. 3.7
end

% Reconstrução com 5 e 10 harmônicas
n_componentes = [5 10 30 100];                 % Número de harmônicas em cada subplot

for k = 1:length(n_componentes)
    x1 = zeros(1,N);                    % Reinicia a reconstrução
    for m = 1:n_componentes(k)
        x1 = x1 + X_mag(m)*cos(2*pi*f(m)*t + X_phase(m));   % Eq. 3.15
    end
    x1 = x1 + a0/2;                     % Adiciona a média

    subplot(length(n_componentes), 1, k);
    plot(t, x1, 'k'); hold on;
    plot(t, x, '--k'); hold off;
    xlabel('Time (s)');
    ylabel('X(t)');
    axis([0 1 -0.2 0.6]);
    text(0.75, 0.5, sprintf('%d components', n_componentes(k)));
end

% A figura mostra a reconstrução de um sinal pela série de Fourier, usando apenas os primeiros 
% harmônicos. 

% O sinal original (linha tracejada)

% É um sinal periódico com período de 1 s, que tem duas partes:

%     - de 0 a 0,5 s, uma rampa que sobe linearmente de 0 até 0,5;
%     - em 0,5 s, uma queda abrupta para zero;
%     - de 0,5 a 1 s, o sinal fica constante em zero.

% Como o período é 1 s, a frequência fundamental é f₀ = 1 Hz, e os harmônicos estão em 1,2,3,4... Hz.

% As aproximações (linha contínua)

% Ao fazer o loop com 5, o sinal é reconstruído somando apenas os 5 primeiros harmônicos (até 5 Hz). 
% e ao fazer o lopp com 10, o sinal é recontruído com 10 harmônicos (até 10 Hz).

% Na rampa, a parte suave do sinal, a aproximação já é boa com 5 componentes e fica quase 
% perfeita com 10. Trechos suaves dependem principalmente das frequências baixas, que já estão 
% incluídas.

% Na descontinuidade (0,5 s), a dificuldade é maior. Uma queda instantânea exige frequências 
% infinitamente altas para ser reproduzida exatamente. Com 5 componentes, a curva desce de forma 
% lenta e arredondada. Com 10, a descida é bem mais íngreme, mas ainda não é vertical. 
% Nos dois casos, a curva cruza 0,5 s exatamente em 0,25, que é a média entre os valores 
% antes (0,5) e depois (0) do salto. 
% 
% É isso que o teorema de Fourier prevê para pontos de descontinuidade.

% É possivel fazer com mais harmonicas?
%
% Sim. O único limite é a frequência de Nyquist: com fs = 500 Hz e fundamental de 1 Hz, 
% você pode usar até 250 harmônicas, porque fs/2 = 250 Hz. 
% Acima disso as senoides não são mais representadas corretamente pelas amostras, e os 
% coeficientes calculados deixam de ter sentido (aliasing).
