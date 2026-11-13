close all;   
clear all;   
clc;         

% Paramêtros do arquivo de dados
% Canal de EMG a ser analisado. 
canal = 1;  % R BICEPS BRAQUIAL
% canal = 11; % L BICEPS BRAQUIAL

% % Envoltório linear: passa-baixa aplicado ao sinal retificado
% % (o texto sugere 3-50 Hz; a FIGURA 4 usa 5 Hz)
% fc_env    = 5;         % frequência de corte (Hz)
% ordem_env = 4;         % ordem do filtro

% % RMS móvel: segundo o SENIAM, para sinais > 50% CVM usar janela de
% % 0,25 a 0,5 s (o típico citado por De Luca é 100-200 ms)
% janela_rms = 0.250;    % duração da janela do RMS (s)

% % Detecção das contrações a partir do envoltório linear
% limiar_rel  = 0.10;    % limiar = 10% do valor máximo do envoltório
% dur_min     = 1.0;     % duração mínima de uma contração (s)
% margem_psd  = 0.5;     % tempo descartado no início/fim de cada contração
%                        % para a análise espectral (s) - evita as rampas

%  ------------------------------------------------------------------------
%  1. LEITURA DO ARQUIVO .HPF
%  ------------------------------------------------------------------------
% arquivo = 'Sub-C005_tempo-0_tarefa-cvim_braço-nt_medida-torque_Plot_and_Store_Rep_1.0.hpf';
% arquivo = 'Sub-C005_tempo-4_tarefa-cvim_braço-nt_medida-emg_Plot_and_Store_Rep_1.0.hpf';

% arquivo = 'Sub-C005_tempo-0_tarefa-cvim_braço-t_medida-torque_Plot_and_Store_Rep_1.0.hpf';
arquivo = 'Sub-C005_tempo-4_tarefa-cvim_braço-t_medida-emg_Plot_and_Store_Rep_1.0.hpf';


hpf = lerHPF(arquivo);            % lê todos os canais do arquivo

% O XML de cada canal traz o nome, a unidade e a frequência de amostragem.
% Extraímos esses campos canal a canal para não depender da ordem do vetor
% hpf.fs (que mistura várias tags com "Frequency"/"Rate" no nome).
xml    = strjoin(hpf.xmlCanais, newline);
nomes  = regexp(xml, '<Name>([^<]*)</Name>', 'tokens');
nomes  = cellfun(@(x) x{1}, nomes, 'UniformOutput', false);
taxas  = regexp(xml, '<PerChannelSampleRate>([^<]*)</PerChannelSampleRate>', 'tokens');
taxas  = cellfun(@(x) str2double(strrep(x{1}, ',', '.')), taxas);
unids  = regexp(xml, '<Unit>([^<]*)</Unit>', 'tokens');
unids  = cellfun(@(x) x{1}, unids, 'UniformOutput', false);

nome_canal = nomes{canal};        % nome do músculo/canal escolhido
fs         = taxas(canal);        % frequência de amostragem (Hz)

fprintf('\nCanal analisado: %d - %s\n', canal, nome_canal);
fprintf('Unidade original: %s | fs = %.2f Hz\n', unids{canal}, fs);

% Sinal de EMG
% emg_bruto = double(hpf.dados{canal});
% Segundo o arquivo o sinal está em volts, por isso converti para milivolts (mV)
emg_bruto = double(hpf.dados{canal}) * 1e3;

% % Remove o nível DC (offset) subtraindo a média do sinal
emg_bruto = emg_bruto - mean(emg_bruto);

% Vetor de tempo (s): cada amostra está separada de 1/fs segundos
N = numel(emg_bruto);
t = (0:N-1)' / fs;
fprintf('Duração do registro: %.1f s (%d amostras)\n', t(end), N);

% Trecho que será utilizado para análise
trecho = [140 160];    % (s)
idx_trecho = t >= trecho(1) & t <= trecho(2);

%  ------------------------------------------------------------------------
%  2. SINAL BRUTO E RETIFICADO E SEU ESPECTRO DE POTÊNCIA (PSD)
%  ------------------------------------------------------------------------
janela_psd = 0.5;      % duração de cada janela de Hamming (s)

% Parâmetros do método de Welch:
%   - janela de Hamming de 'janela_psd' segundos
%   - sobreposição de 50% entre janelas (nwin/2)
%   - nfft = próxima potência de 2 (zero padding) -> resolução = fs/nfft

nwin  = round(janela_psd * fs);           % amostras por janela
nover = round(nwin/2);                    % 50% de sobreposição
nfft  = 2^nextpow2(nwin);                 % número de pontos da FFT
fprintf('PSD (Welch): janela Hamming de %d amostras, nfft = %d, ', nwin, nfft);
fprintf('resolução = %.2f Hz\n', fs/nfft);

% PSD do sinal bruto
[pxx_bruto, f] = pwelch(emg_bruto(idx_trecho), hamming(nwin), nover, nfft, fs);

% Sinal Bruto
subplot(3,1,1);                            % trecho escolhido
plot(t(idx_trecho), emg_bruto(idx_trecho), 'b');
xlabel('Tempo (s)'); ylabel('Amplitude (mV)');
title(['Sinal EMG bruto - ' nome_canal]);
subtitle(sprintf('Trecho de %g a %g s', trecho(1), trecho(2)));
grid on;

% PSD do sinal retificado
emg_retificado = abs(emg_bruto);
subplot(3,1,2);                            % trecho escolhido
plot(t(idx_trecho), emg_retificado(idx_trecho), 'b');
xlabel('Tempo (s)'); ylabel('Amplitude (mV)');
title(['Sinal EMG bruto - ' nome_canal]);
subtitle(sprintf('Trecho de %g a %g s', trecho(1), trecho(2)));
grid on;

subplot(3,1,3);                            % PSD do sinal bruto
plot(f, pxx_bruto, 'r');
xlabel('Frequência (Hz)'); ylabel('PSD (mV^2/Hz)');
title('PSD do sinal bruto');
subtitle(sprintf('Trecho de %g a %g s', trecho(1), trecho(2)));
xlim([0 fs/2]); grid on;

% ------------------------------------------------------------------------
%  3. FILTRAGEM PASSA-BANDA (Butterworth 20-450 Hz, 4a ordem)
%  ------------------------------------------------------------------------
% % Filtro passa-banda (recomendação DELSYS: 20-450 Hz)
% fc_pb    = [20 450];   % frequências de corte inferior e superior (Hz)
% ordem_pb = 4;          % ordem do filtro Butterworth

% % Filtro rejeita-banda para a rede elétrica (60 Hz no Brasil)
% usar_notch = false;    % true para ativar
% f_rede     = 60;       % frequência da rede elétrica (Hz)