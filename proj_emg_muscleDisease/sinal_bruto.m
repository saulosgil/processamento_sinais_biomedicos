%% Sinal bruto de EMG - PhysioNet EMGDB (https://physionet.org/content/emgdb/1.0.0/)
% Sinais de EMG de agulha do músculo tibial anterior:
%   - emg_healthy.txt    : indivíduo saudável
%   - emg_myopathy.txt   : paciente com miopatia
%   - emg_neuropathy.txt : paciente com neuropatia
% Formato dos arquivos: coluna 1 = tempo (s), coluna 2 = amplitude (mV)
% Frequência de amostragem: 4 kHz

clear all; 
close all; 
clc;

%% Arquivos e rótulos
pasta    = fileparts(mfilename('fullpath'));
arquivos = {'emg_healthy.txt', 'emg_myopathy.txt', 'emg_neuropathy.txt'};
titulos  = {'Saudável', 'Miopatia', 'Neuropatia'};
cores    = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880]};

nSinais = numel(arquivos);
emg     = struct('nome', titulos, 't', [], 'x', [], 'fs', []);

%% Leitura dos dados
for k = 1:nSinais
    dados = load(fullfile(pasta, arquivos{k}));   % [tempo, amplitude]
    emg(k).t  = dados(:,1);
    emg(k).x  = dados(:,2);
    emg(k).fs = 1 / mean(diff(emg(k).t));          % ~4000 Hz

    fprintf('%-10s | N = %6d amostras | fs = %.0f Hz | duração = %.2f s\n', ...
        titulos{k}, numel(emg(k).x), emg(k).fs, emg(k).t(end));
end


fs = emg(1).fs;  % frequência de amostragem (Hz)
t = emg(1).t;   % vetor de tempo (s)




% %% Figura 1: sinais brutos completos
% figure('Name', 'EMG - sinal bruto completo', 'Color', 'w');
% for k = 1:nSinais
%     ax(k) = subplot(nSinais, 1, k); %#ok<SAGROW>
%     plot(emg(k).t, emg(k).x, 'Color', cores{k});
%     grid on;
%     title(['EMG bruto - ' titulos{k}]);
%     ylabel('Amplitude (mV)');
%     xlim([emg(k).t(1) emg(k).t(end)]);
% end
% xlabel('Tempo (s)');
% linkaxes(ax, 'y');

%  ------------------------------------------------------------------------
%  2. ESPECTRO DE POTÊNCIA (PSD)
%  ------------------------------------------------------------------------
%% PSD - Densidade Espectral de Potência (método de Welch)
% Janela de Hamming de 1024 amostras (~0,256 s), 50% de sobreposição
% Resolução em frequência: fs/nfft = 4000/1024 ≈ 3,9 Hz
nJanela   = 1024;
nOverlap  = nJanela/2;
nfft      = 1024;

for k = 1:nSinais
    x = emg(k).x - mean(emg(k).x);                 % remove o nível DC
    [emg(k).Pxx, emg(k).f] = pwelch(x, hamming(nJanela), nOverlap, nfft, emg(k).fs);

    % Frequências média e mediana do espectro
    P = emg(k).Pxx;  f = emg(k).f;
    emg(k).fmedia   = sum(f .* P) / sum(P);
    Pacum           = cumsum(P) / sum(P);
    emg(k).fmediana = f(find(Pacum >= 0.5, 1, 'first'));

    fprintf('%-10s | f média = %6.1f Hz | f mediana = %6.1f Hz\n', ...
        titulos{k}, emg(k).fmedia, emg(k).fmediana);
end

%% Figura 3: PSD de cada sinal
figure('Name', 'EMG - PSD (Welch)', 'Color', 'w');
for k = 1:nSinais
    subplot(nSinais, 1, k);
    plot(emg(k).f, 10*log10(emg(k).Pxx), 'Color', cores{k}); hold on;
    xline(emg(k).fmediana, '--k', sprintf('f_{med} = %.0f Hz', emg(k).fmediana));
    grid on;
    title(['PSD - ' titulos{k}]);
    ylabel('PSD (dB/Hz)');
    xlim([0 emg(k).fs/2]);
end
xlabel('Frequência (Hz)');

%% Figura 4: PSDs sobrepostas para comparação
figure('Name', 'EMG - PSD comparativa', 'Color', 'w');
hold on;
for k = 1:nSinais
    plot(emg(k).f, 10*log10(emg(k).Pxx), 'Color', cores{k}, 'LineWidth', 1.2);
end
grid on;
title('PSD - comparação entre os sinais');
xlabel('Frequência (Hz)');
ylabel('PSD (dB/Hz)');
xlim([0 emg(1).fs/2]);
legend(titulos, 'Location', 'northeast');
