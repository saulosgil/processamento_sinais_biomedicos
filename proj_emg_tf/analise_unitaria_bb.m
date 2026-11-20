%% ========================================================================
%  ANÁLISE DE EMG DE UM ARQUIVO .HPF (CVIM)
%  Leitura -> pré-processamento -> PSD -> filtragem -> trecho -> envoltórias -> RMS pico
%  ========================================================================
close all;
clear all;
clc;

%% 1. PARÂMETROS
% ---- Arquivo e canal ----
arquivo = 'Sub-C005_tempo-4_tarefa-cvim_braço-t_medida-emg_Plot_and_Store_Rep_1.0.hpf';

% Escolha o canal desejado: 1 = R BICEPS BRAQUIAL | 21 = R TRICEPS BRAQUIAL
canal   = 1;   

% ---- PSD (método de Welch) ----
janela_psd = 0.5;       % duração de cada janela de Hamming (s)

% ---- Filtragem ----
f_rede   = 60;          % frequência da rede elétrica (Hz)
bw_rede  = 2;           % meia-largura do rejeita-faixa em torno da rede (Hz)
fc_pb    = [20 450];    % cortes do passa-banda (Hz)
ordem_pb = 4;           % ordem final do Butterworth passa-banda

% ---- Análise ----
t_desc   = 3;           % segundos iniciais descartados (transientes)
janela_s = 0.200;       % janela das envoltórias (s)

%% 2. LEITURA DO ARQUIVO HPF
hpf = lerHPF(arquivo);

% O XML de cada canal traz nome e frequência de amostragem. Os campos são
% extraídos canal a canal para não depender da ordem de hpf.fs (que mistura
% várias tags com "Frequency"/"Rate" no nome).
xml   = strjoin(hpf.xmlCanais, newline);
nomes = regexp(xml, '<Name>([^<]*)</Name>', 'tokens');
nomes = cellfun(@(x) x{1}, nomes, 'UniformOutput', false);
taxas = regexp(xml, '<PerChannelSampleRate>([^<]*)</PerChannelSampleRate>', 'tokens');
taxas = cellfun(@(x) str2double(strrep(x{1}, ',', '.')), taxas);

nome_canal = nomes{canal};      % nome do músculo/canal escolhido
fs         = taxas(canal);      % frequência de amostragem (Hz)

%% 3. PRÉ-PROCESSAMENTO
emg_bruto = double(hpf.dados{canal}) * 1e3;   % V -> mV
emg_bruto = emg_bruto - mean(emg_bruto);      % remove o nível DC (offset)

N = numel(emg_bruto);
t = (0:N-1)' / fs;                            % vetor de tempo (s)

% Remove picos isolados: compara cada amostra com a mediana local
% (janela de ~50 ms) e substitui as que passam de 6 desvios (MAD)
jan_out   = round(0.05 * fs);
emg_bruto = filloutliers(emg_bruto, 'linear', 'movmedian', jan_out, ...
                         'ThresholdFactor', 6);
%% 4. FILTRAGEM
% 4.1 Rejeita-faixa em torno de 60 Hz (rede elétrica)
[b_n, a_n]    = butter(2, [f_rede-bw_rede f_rede+bw_rede] / (fs/2), 'stop');
emg_filt_rede = filtfilt(b_n, a_n, emg_bruto);

% 4.2 Passa-banda 20-450 Hz
Wn           = fc_pb / (fs/2);                % normaliza pela frequência de Nyquist (o butter exige valores entre 0 e 1)
[b_pb, a_pb] = butter(ordem_pb/2, Wn, 'bandpass');   % no modo bandpass o butter dobra a ordem -> ordem 4
emg_filt     = filtfilt(b_pb, a_pb, emg_filt_rede);  % filtfilt: fase zero

%% 5. ESPECTRO DE POTÊNCIA (Welch)
%   - janela de Hamming de 'janela_psd' segundos
%   - sobreposição de 50% entre janelas
%   - nfft = próxima potência de 2 (zero padding) -> resolução = fs/nfft
nwin  = round(janela_psd * fs);
nover = round(nwin/2);
nfft  = 2^nextpow2(nwin);

[pxx_bruto, f]  = pwelch(emg_bruto,     hamming(nwin), nover, nfft, fs);
pxx_filt_rede   = pwelch(emg_filt_rede, hamming(nwin), nover, nfft, fs);
pxx_filt        = pwelch(emg_filt,      hamming(nwin), nover, nfft, fs);

%% 6. TRECHO DE ANÁLISE (exclui os primeiros t_desc segundos)
trecho     = [min(t) + t_desc, max(t)];             % (s)
idx_trecho = t >= trecho(1) & t <= trecho(2);       % índice lógico

t_trecho        = t(idx_trecho);
emg_filt_trecho = emg_filt(idx_trecho);

%% 7. RETIFICAÇÃO, ENVOLTÓRIAS E RMS PICO
janela = round(janela_s * fs);                      % janela em amostras

emg_retificado = abs(emg_filt_trecho);
env_movmean    = movmean(emg_retificado, janela);           % envoltória por média móvel
env_rms        = sqrt(movmean(emg_filt_trecho.^2, janela)); % envoltória RMS móvel centrada

[rms_pico, i_max] = max(env_rms);                   % maior RMS dentro do trecho

fprintf('Arquivo : %s\n', arquivo);
fprintf('Canal   : %d - %s | fs = %.2f Hz\n', canal, nome_canal, fs);
fprintf('RMS pico (janela %d ms): %.4f mV em t = %.2f s\n', ...
        round(janela_s*1000), rms_pico, t_trecho(i_max));

%% 8. FIGURAS
ms = round(janela_s*1000);

% 8.1 Sinal bruto e seu PSD
figure('Name', 'Sinal bruto');
subplot(2,1,1);
plot(t, emg_bruto, 'b');
xlabel('Tempo (s)'); ylabel('Amplitude (mV)');
title(['Sinal EMG bruto - ' nome_canal]); grid on;

subplot(2,1,2);
plot(f, pxx_bruto, 'r', 'LineWidth', 1.5);
xlabel('Frequência (Hz)'); ylabel('PSD (mV^2/Hz)');
title('PSD do sinal bruto'); grid on;

% 8.2 PSD em cada etapa da filtragem
figure('Name', 'PSD por etapa');
ax1 = subplot(2,3,1);
plot(f, pxx_bruto, 'r', 'LineWidth', 1.5);
title('Bruto'); xlabel('Frequência (Hz)'); ylabel('PSD (mV^2/Hz)'); grid on;

ax2 = subplot(2,3,2);
plot(f, pxx_filt_rede, 'b', 'LineWidth', 1.5);
title('Notch 60 Hz'); xlabel('Frequência (Hz)'); ylabel('PSD (mV^2/Hz)'); grid on;

ax3 = subplot(2,3,3);
plot(f, pxx_filt, 'g:', 'LineWidth', 1.5);
title('Notch + passa-banda'); xlabel('Frequência (Hz)'); ylabel('PSD (mV^2/Hz)'); grid on;

linkaxes([ax1 ax2 ax3], 'xy');      % mesma escala nos três
xlim(ax1, [0 fs/2]);

subplot(2,3,[4 5 6]);
plot(f, pxx_bruto,     'r',  'LineWidth', 1.5); hold on;
plot(f, pxx_filt_rede, 'b',  'LineWidth', 1.5);
plot(f, pxx_filt,      'g:', 'LineWidth', 1.5); hold off;
xline(f_rede, ':k', '60 Hz');
xlabel('Frequência (Hz)'); ylabel('PSD (mV^2/Hz)');
legend('Bruto', 'Notch 60 Hz', 'Notch + passa-banda 20–450 Hz');
xlim([0 fs/2]); grid on;

% 8.3 Sinal bruto x filtrado no tempo
figure('Name', 'Bruto x filtrado');
plot(t, emg_bruto, 'b', 'LineWidth', 1); hold on;
plot(t, emg_filt,  'r', 'LineWidth', 1); hold off;
xline(trecho(1), '--k', 'Início do trecho');
xlabel('Tempo (s)'); ylabel('Amplitude (mV)');
legend('Sinal bruto', 'Sinal filtrado (Notch + passa-banda 20–450 Hz)');
grid on;

% 8.4 Trecho analisado: retificado, envoltórias e RMS pico
figure('Name', 'Envoltórias');
plot(t_trecho, emg_retificado, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5); hold on;
plot(t_trecho, env_movmean, 'b', 'LineWidth', 2);
plot(t_trecho, env_rms,     'r', 'LineWidth', 2);
plot(t_trecho(i_max), rms_pico, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 6); hold off; % marca o RMS pico
xlabel('Tempo (s)'); ylabel('Amplitude (mV)');
title(nome_canal);
legend('EMG filtrado retificado', ...
       sprintf('Envoltória média móvel (%d ms)', ms), ...
       sprintf('Envoltória RMS (%d ms)', ms), ...
       sprintf('RMS pico = %.4f mV', rms_pico));
grid on;
%%
% Códigos antes de serem organizados pela IA comentados
% clear all;
% clc;

% %% Leitura do arquivo HPF
% % Canal de EMG a ser analisado.
% canal = 1;  % R BICEPS BRAQUIAL
% % canal = 21;  % R TRICEPS BRAQUIAL

% % caminho do arquivo de dados
% arquivo = 'Sub-C001_tempo-0_tarefa-cvim_braço-t_medida-emg_Plot_and_Store_Rep_1.0.hpf';

% % Lê o arquivo HPF
% hpf = lerHPF(arquivo);

% % O XML de cada canal traz o nome, a unidade e a frequência de amostragem.
% % Extrai os campos canal a canal para não depender da ordem do vetor
% % hpf.fs (que mistura várias tags com "Frequency"/"Rate" no nome).
% xml    = strjoin(hpf.xmlCanais, newline);
% nomes  = regexp(xml, '<Name>([^<]*)</Name>', 'tokens');
% nomes  = cellfun(@(x) x{1}, nomes, 'UniformOutput', false);
% taxas  = regexp(xml, '<PerChannelSampleRate>([^<]*)</PerChannelSampleRate>', 'tokens');
% taxas  = cellfun(@(x) str2double(strrep(x{1}, ',', '.')), taxas);
% % unids  = regexp(xml, '<Unit>([^<]*)</Unit>', 'tokens');
% % unids  = cellfun(@(x) x{1}, unids, 'UniformOutput', false);

% nome_canal = nomes{canal};        % nome do músculo/canal escolhido
% fs         = taxas(canal);        % frequência de amostragem (Hz)

% % % fprintf('\nCanal analisado: %d - %s\n', canal, nome_canal); % comentado para não poluir o log do script
% % fprintf('fs = %.2f Hz\n', fs);                                % comentado para não poluir o log do script

% %% Pré-processamento do sinal de EMG
% emg_bruto = double(hpf.dados{canal}) * 1e3; % Converte para mV
% emg_bruto = emg_bruto - mean(emg_bruto);    % Remove o nível DC (offset) subtraindo a média do sinal

% % Vetor de tempo (s): cada amostra está separada de 1/fs segundos
% N = numel(emg_bruto);
% t = (0:N-1)' / fs;

% %% Espectro de potência (PSD) usando o método de Welch
% % Parâmetros do método de Welch:
%     %   - janela de Hamming de 'janela_psd' segundos
%     %   - sobreposição de 50% entre janelas (nwin/2)
%     %   - nfft = próxima potência de 2 (zero padding) -> resolução = fs/nfft
% janela_psd = 0.5;                         % duração de cada janela de Hamming (s)
% nwin  = round(janela_psd * fs);           % amostras por janela
% nover = round(nwin/2);                    % 50% de sobreposição
% nfft  = 2^nextpow2(nwin);                 % número de pontos da FFT

% % PSD do sinal bruto
% [pxx_bruto, f] = pwelch(emg_bruto, hamming(nwin), nover, nfft, fs);

% %% Gráficos do sinal bruto e do PSD
% % Sinal bruto
% figure;
% subplot(2,1,1);
% plot(t, emg_bruto);
% xlabel('Tempo (s)');
% ylabel('Amplitude (mV)');
% grid on;

% % PSD do sinal bruto
% subplot(2,1,2);
% plot(f, pxx_bruto, 'r', 'LineWidth', 1.5);
% xlabel('Frequência (Hz)');
% ylabel('PSD (dB/Hz)');
% grid on;

% %% Filtragem
% % Rede eletrica  - notch filter 60Hz
% f_rede = 60;                                                                    % Frequência da rede elétrica
% [b_n, a_n] = butter(2, [f_rede-2 f_rede+2] / (fs/2), 'stop');
% emg_filt_rede = filtfilt(b_n, a_n, emg_bruto);

% % PSD do sinal filtrado
% [pxx_filt_rede, f] = pwelch(emg_filt_rede, hamming(nwin), nover, nfft, fs);

% % Filtro passa-banda (recomendação: 20-450 Hz)
% fc_pb    = [20 450];                                                            % frequências de corte inferior e superior (Hz)
% ordem_pb = 4;                                                                   % ordem do filtro Butterworth
% Wn = fc_pb / (fs/2);                                                            % normaliza as frequências de corte pela frequência de Nyquist (fs/2), pois o butter exige valores entre 0 e 1
% [b_pb, a_pb] = butter(ordem_pb/2, Wn, 'bandpass');

% % filtfilt - aplica o filtro
% emg_filt = filtfilt(b_pb, a_pb, emg_filt_rede);

% % PSD do sinal filtrado
% [pxx_filt, f] = pwelch(emg_filt, hamming(nwin), nover, nfft, fs);

% figure;
% subplot(2,3,1);
% plot(f, pxx_bruto, 'r', 'LineWidth', 1.5);
% xlabel('Frequência (Hz)');
% ylabel('PSD (dB/Hz)');
% grid on;

% subplot(2,3,2);
% plot(f, pxx_filt_rede, 'b', 'LineWidth', 1.5); hold on;
% xlabel('Frequência (Hz)');
% ylabel('PSD (dB/Hz)');
% grid on;

% subplot(2,3,3);
% plot(f, pxx_filt, 'g:', 'LineWidth', 1.5);
% xlabel('Frequência (Hz)');
% ylabel('PSD (dB/Hz)');
% grid on;

% subplot(2,3,[4 5 6]);
% plot(f, pxx_bruto, 'r', 'LineWidth', 1.5); hold on;
% plot(f, pxx_filt_rede, 'b', 'LineWidth', 1.5); hold on;
% plot(f, pxx_filt, 'g:', 'LineWidth', 1.5);
% xlabel('Frequência (Hz)');
% ylabel('PSD (dB/Hz)');
% legend('Bruto', 'Rede - 60 Hz (Notch)', 'Notch + passa-banda 20–450 Hz');
% grid on;

% figure;
% plot(t, emg_bruto, 'b', 'LineWidth', 1.5); hold on;
% plot(t, emg_filt, 'r', 'LineWidth', 1.5);
% xlabel('Tempo (s)');
% ylabel('Amplitude (mV)');
% legend('Sinal bruto', 'Sinal filtrado (Notch + passa-banda 20–450 Hz)');
% grid on;

% %% Trecho que será utilizado para análise - exluindo os primeiros 3 segundos (transientes)
% inicio_trecho = min(t) + 3;                      % (s)
% fim_trecho    = max(t);                          % (s)
% trecho = [inicio_trecho, fim_trecho];            % (s)

% idx_trecho = t >= trecho(1) & t <= trecho(2);    % índice lógico: true nas amostras dentro do trecho

% t_trecho        = t(idx_trecho);
% emg_filt_trecho = emg_filt(idx_trecho);

% figure;
% plot(t_trecho, emg_filt_trecho, 'r', 'LineWidth', 1.5);
% xlabel('Tempo (s)');
% ylabel('Amplitude (mV)');
% title('Sinal filtrado (Notch + passa-banda 20–450 Hz)');
% legend('Sinal filtrado (Notch + passa-banda 20–450 Hz)');
% grid on;

% %% Sinal retificado e sinal com envoltória (média móvel e RMS)
% janela_s = 0.200;                                           % janela de 200 ms
% janela   = round(janela_s * fs);                            % em amostras

% emg_retificado = abs(emg_filt_trecho);
% emg_movmean = movmean(emg_retificado, janela);              % envoltória do sinal retificado (média móvel)
% env_rms = sqrt(movmean(emg_filt_trecho.^2, janela));        % RMS móvel centrado

% figure;
% plot(t_trecho, emg_retificado, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5); hold on;
% plot(t_trecho, emg_movmean, 'b', 'LineWidth', 2); hold on;
% plot(t_trecho, env_rms, 'r', 'LineWidth', 2); hold off;
% xlabel('Tempo (s)');
% ylabel('Amplitude (mV)');
% legend('EMG filtrado', sprintf('Envoltória média móvel (%d ms)', janela_s*1000), sprintf('Envoltória RMS (%d ms)', janela_s*1000));
% grid on;

% % Encontra o maior RMS dentro do trecho
% [rms_cvim, i_max] = max(env_rms);                           % maior RMS dentro do trecho
% fprintf('RMS pico (janela 200 ms): %.4f mV em t = %.2f s\n', rms_cvim, t_trecho(i_max));