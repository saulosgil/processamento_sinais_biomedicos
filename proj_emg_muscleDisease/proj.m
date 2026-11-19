clear all; 
close all; 
clc;

%% Indivíduo saudável
arquivo = 'emg_healthy.txt';    % sinal de EMG de um indivíduo saudável
dados   = load(arquivo);        % [tempo, amplitude]
t = dados(:,1);                 % vetor de tempo (s)
x = dados(:,2);                 % sinal de EMG (mV)
fs = 1 / mean(diff(t));         % frequência de amostragem (Hz)

subplot(3, 2, 1);
plot(t, x, 'b');
grid on;
title('Sinal de EMG bruto - Indivíduo saudável');
xlabel('Tempo (s)');
ylabel('Amplitude (mV)');

% Espectro de potencia do sinal de EMG bruto
% Parâmetros do método de Welch:
%   - janela de Hamming de 'janela_psd' segundos
%   - sobreposição de 50% entre janelas (nwin/2)
%   - nfft = próxima potência de 2 (zero padding) -> resolução = fs/nfft

janela_psd = 0.5;                         % duração de cada janela de Hamming (s)
nwin  = round(janela_psd * fs);           % amostras por janela
nover = round(nwin/2);                    % 50% de sobreposição
nfft  = 2^nextpow2(nwin);                 % número de pontos da FFT
fprintf('PSD (Welch): janela Hamming de %d amostras, nfft = %d, ', nwin, nfft);
fprintf('resolução = %.2f Hz\n', fs/nfft);

% PSD do sinal bruto
[pxx_bruto, f] = pwelch(x, hamming(nwin), nover, nfft, fs);

subplot(3, 2, 2);
plot(f, pxx_bruto, 'r');
xlabel('Frequência (Hz)'); ylabel('PSD (mV^2/Hz)');
title('PSD do sinal bruto');
xlim([0 400]); 
grid on;

%% Indivíduo miopatia
arquivo = 'emg_myopathy.txt';   % sinal de EMG de um paciente com miopatia
dados   = load(arquivo);        % [tempo, amplitude]
t = dados(:,1);                 % vetor de tempo (s)
x = dados(:,2);                 % sinal de EMG (mV)
fs = 1 / mean(diff(t));         % frequência de amostragem (Hz)

subplot(3, 2, 3);
plot(t, x, 'g');
grid on;
title('Sinal de EMG bruto - Indivíduo com miopatia');
xlabel('Tempo (s)');
ylabel('Amplitude (mV)');

% Espectro de potencia do sinal de EMG bruto
% Parâmetros do método de Welch:
%   - janela de Hamming de 'janela_psd' segundos
%   - sobreposição de 50% entre janelas (nwin/2)
%   - nfft = próxima potência de 2 (zero padding) -> resolução = fs/nfft

janela_psd = 0.5;                         % duração de cada janela de Hamming (s)
nwin  = round(janela_psd * fs);           % amostras por janela
nover = round(nwin/2);                    % 50% de sobreposição
nfft  = 2^nextpow2(nwin);                 % número de pontos da FFT
fprintf('PSD (Welch): janela Hamming de %d amostras, nfft = %d, ', nwin, nfft);
fprintf('resolução = %.2f Hz\n', fs/nfft);

% PSD do sinal bruto
[pxx_bruto, f] = pwelch(x, hamming(nwin), nover, nfft, fs);

subplot(3, 2, 4);
plot(f, pxx_bruto, 'r');
xlabel('Frequência (Hz)'); ylabel('PSD (mV^2/Hz)');
title('PSD do sinal bruto');
xlim([0 400]); 
grid on;

%% Indivíduo neuropatia
arquivo = 'emg_neuropathy.txt'; % sinal de EMG de um paciente com neuropatia
dados   = load(arquivo);        % [tempo, amplitude]
t = dados(:,1);                 % vetor de tempo (s)
x = dados(:,2);                 % sinal de EMG (mV)
fs = 1 / mean(diff(t));         % frequência de amostragem (Hz)

subplot(3, 2, 5);
plot(t, x, 'r');
grid on;
title('Sinal de EMG bruto - Indivíduo com neuropatia');
xlabel('Tempo (s)');
ylabel('Amplitude (mV)');

% Espectro de potencia do sinal de EMG bruto
% Parâmetros do método de Welch:
%   - janela de Hamming de 'janela_psd' segundos
%   - sobreposição de 50% entre janelas (nwin/2)
%   - nfft = próxima potência de 2 (zero padding) -> resolução = fs/nfft

janela_psd = 0.5;                         % duração de cada janela de Hamming (s)
nwin  = round(janela_psd * fs);           % amostras por janela
nover = round(nwin/2);                    % 50% de sobreposição
nfft  = 2^nextpow2(nwin);                 % número de pontos da FFT
fprintf('PSD (Welch): janela Hamming de %d amostras, nfft = %d, ', nwin, nfft);
fprintf('resolução = %.2f Hz\n', fs/nfft);

% PSD do sinal bruto
[pxx_bruto, f] = pwelch(x, hamming(nwin), nover, nfft, fs);

subplot(3, 2, 6);
plot(f, pxx_bruto, 'r');
xlabel('Frequência (Hz)'); ylabel('PSD (mV^2/Hz)');
title('PSD do sinal bruto');
xlim([0 400]); 
grid on;