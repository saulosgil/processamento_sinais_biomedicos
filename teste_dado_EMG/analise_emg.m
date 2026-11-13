%% ========================================================================
%  ANÁLISE DO SINAL EMG - CONTRAÇÃO VOLUNTÁRIA ISOMÉTRICA MÁXIMA (CVIM)
%  ------------------------------------------------------------------------
%  Baseado em: Marchetti & Duarte (2006) "Instrumentação em Eletromiografia"
%
%  Etapas (seções do texto entre parênteses):
%    1. Leitura do arquivo .hpf e escolha do canal de EMG
%    2. Sinal bruto e PSD do sinal bruto                     (1.1.4, 1.2.3)
%    3. Filtragem passa-banda Butterworth 20-450 Hz, 4a ordem,
%       com fase zero (filtfilt)                             (1.1.4)
%    4. Retificação de onda completa                         (1.2.1.1)
%    5. Envoltório linear (passa-baixa 5 Hz)                 (1.2.1.2)
%    6. RMS móvel                                            (1.2.1.3)
%    7. EMG integrado (iEMG)                                 (1.2.1.4)
%    8. Detecção automática das contrações
%    9. Normalização pela CVIM (% CVIM)                      (1.2.2)
%   10. Análise no domínio da frequência: PSD, frequência média
%       (MNF) e frequência mediana (MDF) de cada contração   (1.2.3)
%
%  Requer: Signal Processing Toolbox (butter, filtfilt, pwelch)
%          e a função lerHPF.m na mesma pasta.
%  ========================================================================

close all;   % fecha todas as figuras abertas
clear all;   % limpa todas as variáveis do workspace
clc;         % limpa a janela de comandos

%% ------------------------------------------------------------------------
%  PARÂMETROS DA ANÁLISE (altere aqui)
%  ------------------------------------------------------------------------
arquivo = 'Sub-C005_tempo-0_tarefa-cvim_braço-nt_medida-torque_Plot_and_Store_Rep_1.0.hpf';

% Canal de EMG a ser analisado. Neste arquivo os canais de EMG são:
%   1 = R BICEPS BRACHII (EMG 1)     11 = L BICEPS BRACHII (EMG 2)
%  21 = R TRICEPS BRACHII (EMG 3)    31 = L TRICEPS BRACHII (EMG 4)
% Os demais canais são acelerômetro, giroscópio e magnetômetro (IMU).
% O canal 11 é o que apresenta as contrações máximas (3 tentativas).
canal = 11;

% Filtro passa-banda (recomendação DELSYS: 20-450 Hz)
fc_pb    = [20 450];   % frequências de corte inferior e superior (Hz)
ordem_pb = 4;          % ordem do filtro Butterworth

% Filtro rejeita-banda opcional para a rede elétrica (60 Hz no Brasil)
usar_notch = false;    % true para ativar
f_rede     = 60;       % frequência da rede elétrica (Hz)

% Envoltório linear: passa-baixa aplicado ao sinal retificado
% (o texto sugere 3-50 Hz; a FIGURA 4 usa 5 Hz)
fc_env    = 5;         % frequência de corte (Hz)
ordem_env = 4;         % ordem do filtro

% RMS móvel: segundo o SENIAM, para sinais > 50% CVM usar janela de
% 0,25 a 0,5 s (o típico citado por De Luca é 100-200 ms)
janela_rms = 0.250;    % duração da janela do RMS (s)

% Detecção das contrações a partir do envoltório linear
limiar_rel  = 0.10;    % limiar = 10% do valor máximo do envoltório
dur_min     = 1.0;     % duração mínima de uma contração (s)
margem_psd  = 0.5;     % tempo descartado no início/fim de cada contração
                       % para a análise espectral (s) - evita as rampas

% Estimativa espectral (método de Welch)
janela_psd = 0.5;      % duração de cada janela de Hamming (s)

% Trecho para visualização ampliada (mesmo intervalo usado no test.m)
trecho = [140 160];    % (s)

%% ------------------------------------------------------------------------
%  1. LEITURA DO ARQUIVO .HPF
%  ------------------------------------------------------------------------
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

% Verificação do teorema de Nyquist (seção 1.1.1): fs deve ser pelo menos
% o dobro da maior frequência de interesse (corte superior do filtro)
if fs < 2*fc_pb(2)
    error('fs = %.1f Hz é insuficiente para o corte de %d Hz.', fs, fc_pb(2));
end

% Sinal de EMG: o arquivo está em volts; convertemos para milivolts (mV)
emg_bruto = double(hpf.dados{canal}) * 1e3;

% Remove o nível DC (offset) subtraindo a média do sinal
emg_bruto = emg_bruto - mean(emg_bruto);

% Vetor de tempo (s): cada amostra está separada de 1/fs segundos
N = numel(emg_bruto);
t = (0:N-1)' / fs;
fprintf('Duração do registro: %.1f s (%d amostras)\n', t(end), N);

% Índices do trecho de visualização ampliada
idx_trecho = t >= trecho(1) & t <= trecho(2);

%% ------------------------------------------------------------------------
%  2. SINAL BRUTO E SEU ESPECTRO DE POTÊNCIA (PSD)
%  ------------------------------------------------------------------------
% Parâmetros do método de Welch:
%   - janela de Hamming de 'janela_psd' segundos
%   - sobreposição de 50% entre janelas
%   - nfft = próxima potência de 2 (zero padding) -> resolução = fs/nfft
nwin  = round(janela_psd * fs);           % amostras por janela
nover = round(nwin/2);                    % 50% de sobreposição
nfft  = 2^nextpow2(nwin);                 % número de pontos da FFT
fprintf('PSD (Welch): janela Hamming de %d amostras, nfft = %d, ', nwin, nfft);
fprintf('resolução = %.2f Hz\n', fs/nfft);

% PSD do sinal bruto (usamos apenas o trecho escolhido, com contração)
[pxx_bruto, f] = pwelch(emg_bruto(idx_trecho), hamming(nwin), nover, nfft, fs);

figure('Name', 'Sinal bruto', 'Color', 'w');
subplot(2,2,[1 2]);                        % sinal bruto completo
plot(t, emg_bruto, 'b');
xlabel('Tempo (s)'); ylabel('Amplitude (mV)');
title(['Sinal EMG bruto - ' nome_canal]);
xlim([0 t(end)]); grid on;

subplot(2,2,3);                            % zoom no trecho escolhido
plot(t(idx_trecho), emg_bruto(idx_trecho), 'b');
xlabel('Tempo (s)'); ylabel('Amplitude (mV)');
title(sprintf('Trecho de %g a %g s', trecho(1), trecho(2)));
grid on;

subplot(2,2,4);                            % PSD do sinal bruto
plot(f, pxx_bruto, 'r');
xlabel('Frequência (Hz)'); ylabel('PSD (mV^2/Hz)');
title('PSD do sinal bruto (trecho)');
xlim([0 fs/2]); grid on;

%% ------------------------------------------------------------------------
%  3. FILTRAGEM PASSA-BANDA (Butterworth 20-450 Hz, 4a ordem)
%  ------------------------------------------------------------------------
% O Butterworth tem resposta plana na banda de passagem (sem ripple),
% sendo o filtro recomendado para EMG (seção 1.1.4).
% A função butter recebe as frequências de corte normalizadas pela
% frequência de Nyquist (fs/2). Para um passa-banda, butter(n, ...)
% gera um filtro de ordem 2n; por isso usamos ordem_pb/2.
Wn = fc_pb / (fs/2);
[b_pb, a_pb] = butter(ordem_pb/2, Wn, 'bandpass');

% filtfilt aplica o filtro para frente e depois de trás para frente,
% eliminando o atraso de fase (seção 1.1.4, pág. 14).
emg_filt = filtfilt(b_pb, a_pb, emg_bruto);

% Filtro rejeita-banda opcional para remover a interferência da rede
if usar_notch
    [b_n, a_n] = butter(2, [f_rede-2 f_rede+2] / (fs/2), 'stop');
    emg_filt   = filtfilt(b_n, a_n, emg_filt);
end

% PSD do sinal filtrado no mesmo trecho (para comparar com o bruto)
pxx_filt = pwelch(emg_filt(idx_trecho), hamming(nwin), nover, nfft, fs);

% Resposta em frequência do filtro projetado (magnitude em dB)
[H, fH] = freqz(b_pb, a_pb, 2048, fs);

figure('Name', 'Filtragem', 'Color', 'w');
subplot(2,2,1);                            % sinal filtrado (trecho)
plot(t(idx_trecho), emg_filt(idx_trecho), 'b');
xlabel('Tempo (s)'); ylabel('Amplitude (mV)');
title(sprintf('Sinal filtrado Butterworth (%d-%d Hz)', fc_pb));
grid on;

subplot(2,2,2);                            % PSD: bruto x filtrado
plot(f, pxx_bruto, 'Color', [0.7 0.7 0.7]); hold on;
plot(f, pxx_filt, 'r');
xlabel('Frequência (Hz)'); ylabel('PSD (mV^2/Hz)');
title('PSD: bruto x filtrado');
legend('Bruto', 'Filtrado'); xlim([0 fs/2]); grid on;

subplot(2,2,[3 4]);                        % resposta do filtro
plot(fH, 20*log10(abs(H)), 'k');
xlabel('Frequência (Hz)'); ylabel('Magnitude (dB)');
title('Resposta em frequência do filtro passa-banda');
ylim([-80 5]); xlim([0 fs/2]); grid on;

%% ------------------------------------------------------------------------
%  4. RETIFICAÇÃO DE ONDA COMPLETA (full-wave)
%  ------------------------------------------------------------------------
% Toma o valor absoluto do sinal: as fases negativas são "rebatidas"
% para o lado positivo, preservando a energia do sinal (seção 1.2.1.1).
emg_ret = abs(emg_filt);

%% ------------------------------------------------------------------------
%  5. ENVOLTÓRIO LINEAR
%  ------------------------------------------------------------------------
% Filtro passa-baixa aplicado ao sinal retificado, suavizando as
% flutuações rápidas e deixando apenas a "envoltória" da amplitude
% (seção 1.2.1.2). Também usamos filtfilt para não atrasar o sinal.
[b_env, a_env] = butter(ordem_env, fc_env / (fs/2), 'low');
emg_env = filtfilt(b_env, a_env, emg_ret);

% O filtro pode gerar pequenos valores negativos nas transições;
% como amplitude não pode ser negativa, eles são zerados.
emg_env(emg_env < 0) = 0;

%% ------------------------------------------------------------------------
%  6. RMS MÓVEL
%  ------------------------------------------------------------------------
% RMS = sqrt( (1/N) * soma(x_i^2) ), calculado em uma janela que desliza
% amostra a amostra (janelas sobrepostas) ao longo do sinal filtrado.
% Não é necessário retificar, pois o sinal é elevado ao quadrado.
n_rms   = round(janela_rms * fs);              % amostras na janela
emg_rms = sqrt(movmean(emg_filt.^2, n_rms));   % média móvel de x^2 -> raiz

%% ------------------------------------------------------------------------
%  7. EMG INTEGRADO (iEMG)
%  ------------------------------------------------------------------------
% Área acumulada sob a curva do sinal retificado (integral numérica pela
% regra dos trapézios). A unidade resultante é mV.s (seção 1.2.1.4).
emg_int = cumtrapz(t, emg_ret);

%% --- Figuras no estilo das FIGURAS 4 e 5 do texto -----------------------
figure('Name', 'Domínio do tempo', 'Color', 'w');
ax(1) = subplot(5,1,1);
plot(t, emg_filt, 'b'); ylabel('mV');
title('Sinal EMG filtrado (20-450 Hz)'); grid on;

ax(2) = subplot(5,1,2);
plot(t, emg_ret, 'b'); ylabel('mV');
title('Sinal EMG retificado (full-wave)'); grid on;

ax(3) = subplot(5,1,3);
plot(t, emg_env, 'b', 'LineWidth', 1.5); ylabel('mV');
title(sprintf('Envoltório linear (passa-baixa %g Hz)', fc_env)); grid on;

ax(4) = subplot(5,1,4);
plot(t, emg_rms, 'b', 'LineWidth', 1.5); ylabel('mV');
title(sprintf('RMS móvel (janela de %g ms)', janela_rms*1e3)); grid on;

ax(5) = subplot(5,1,5);
plot(t, emg_int, 'b', 'LineWidth', 1.5); ylabel('mV.s');
title('EMG integrado (iEMG)'); xlabel('Tempo (s)'); grid on;

linkaxes(ax, 'x');        % zoom em um gráfico aplica-se a todos
xlim(ax(1), [0 t(end)]);

% Mesma figura, mas apenas no trecho ampliado (facilita a visualização)
figure('Name', 'Domínio do tempo - trecho', 'Color', 'w');
subplot(3,1,1);
plot(t(idx_trecho), emg_filt(idx_trecho), 'b'); ylabel('mV');
title('Filtrado'); grid on;

subplot(3,1,2);
plot(t(idx_trecho), emg_ret(idx_trecho), 'b'); hold on;
plot(t(idx_trecho), emg_env(idx_trecho), 'r', 'LineWidth', 2);
ylabel('mV'); title('Retificado (azul) e envoltório linear (vermelho)'); grid on;

subplot(3,1,3);
plot(t(idx_trecho), emg_rms(idx_trecho), 'k', 'LineWidth', 1.5);
ylabel('mV'); xlabel('Tempo (s)'); title('RMS móvel'); grid on;

%% ------------------------------------------------------------------------
%  8. DETECÇÃO AUTOMÁTICA DAS CONTRAÇÕES
%  ------------------------------------------------------------------------
% Uma amostra é considerada "ativa" quando o envoltório linear supera
% uma fração (limiar_rel) do seu valor máximo.
limiar = limiar_rel * max(emg_env);
ativo  = emg_env > limiar;

% Localiza inícios (0 -> 1) e fins (1 -> 0) dos trechos ativos.
% Os zeros extras garantem que trechos encostados nas bordas sejam vistos.
d      = diff([0; ativo; 0]);
inicio = find(d == 1);
fim    = find(d == -1) - 1;

% Descarta trechos curtos demais (picos/artefatos isolados)
validos = (fim - inicio + 1) / fs >= dur_min;
inicio  = inicio(validos);
fim     = fim(validos);
nC      = numel(inicio);
fprintf('\n%d contração(ões) detectada(s) (limiar = %.3f mV)\n', nC, limiar);

if nC == 0
    error('Nenhuma contração detectada. Reduza limiar_rel ou dur_min.');
end

%% ------------------------------------------------------------------------
%  9 e 10. MEDIDAS POR CONTRAÇÃO, NORMALIZAÇÃO E ANÁLISE ESPECTRAL
%  ------------------------------------------------------------------------
% Vetores para guardar os resultados de cada contração
dur      = zeros(nC,1);   % duração (s)
rms_pico = zeros(nC,1);   % maior valor do RMS móvel (mV)
rms_med  = zeros(nC,1);   % RMS do trecho central (mV)
iemg     = zeros(nC,1);   % iEMG da contração (mV.s)
mnf      = zeros(nC,1);   % frequência média (Hz)
mdf      = zeros(nC,1);   % frequência mediana (Hz)
pxx_c    = zeros(numel(f), nC);   % PSD de cada contração

n_marg = round(margem_psd * fs);  % amostras descartadas em cada borda

for k = 1:nC
    % Índices da contração completa
    ic = inicio(k):fim(k);

    % Trecho central (sem as rampas de subida e descida), usado para o RMS
    % médio e para o espectro, que exigem um sinal mais estacionário
    i1 = inicio(k) + n_marg;
    i2 = fim(k)    - n_marg;
    if i2 - i1 < nwin            % se a contração for curta, usa tudo
        i1 = inicio(k); i2 = fim(k);
    end
    ip = i1:i2;

    dur(k)      = (fim(k) - inicio(k) + 1) / fs;
    rms_pico(k) = max(emg_rms(ic));
    rms_med(k)  = sqrt(mean(emg_filt(ip).^2));
    iemg(k)     = trapz(t(ic), emg_ret(ic));

    % PSD da contração pelo método de Welch (janela de Hamming)
    pxx_c(:,k) = pwelch(emg_filt(ip), hamming(nwin), nover, nfft, fs);

    % Considera apenas a banda útil do EMG (a do filtro passa-banda)
    banda = f >= fc_pb(1) & f <= fc_pb(2);
    fb = f(banda);
    pb = pxx_c(banda, k);

    % Frequência média:  MNF = soma(f * P(f)) / soma(P(f))
    mnf(k) = sum(fb .* pb) / sum(pb);

    % Frequência mediana: frequência que divide a potência total em
    % duas metades iguais (potência acumulada = 50%)
    pot_acum = cumsum(pb) / sum(pb);
    mdf(k)   = fb(find(pot_acum >= 0.5, 1, 'first'));
end

% Normalização pela CVIM (seção 1.2.2, método 1): o maior RMS encontrado
% entre as tentativas de contração máxima corresponde a 100%.
[cvim, k_cvim] = max(rms_pico);
emg_norm = 100 * emg_rms / cvim;                 % sinal em % da CVIM
rms_pico_pct = 100 * rms_pico / cvim;            % pico de cada tentativa
rms_med_pct  = 100 * rms_med  / cvim;            % média de cada tentativa

% Tabela de resultados na janela de comandos
resultados = table((1:nC)', t(inicio), t(fim), dur, rms_pico, rms_pico_pct, ...
    rms_med, rms_med_pct, iemg, mnf, mdf, ...
    'VariableNames', {'Contracao', 'Inicio_s', 'Fim_s', 'Duracao_s', ...
    'RMSpico_mV', 'RMSpico_pctCVIM', 'RMSmedio_mV', 'RMSmedio_pctCVIM', ...
    'iEMG_mVs', 'MNF_Hz', 'MDF_Hz'});
fprintf('\nValor de referência (CVIM) = %.3f mV (contração %d)\n\n', cvim, k_cvim);
disp(resultados);

%% --- Figura: detecção das contrações e normalização (FIGURAS 6 e 7) -----
figure('Name', 'Contrações e normalização', 'Color', 'w');
subplot(2,1,1);
plot(t, emg_env, 'b', 'LineWidth', 1.2); hold on;
yline(limiar, 'r--', 'Limiar');
for k = 1:nC
    % Sombreia cada contração detectada
    patch([t(inicio(k)) t(fim(k)) t(fim(k)) t(inicio(k))], ...
          [0 0 1 1] * max(emg_env), 'g', ...
          'FaceAlpha', 0.15, 'EdgeColor', 'none');
    text(mean(t([inicio(k) fim(k)])), 1.05*max(emg_env), sprintf('C%d', k), ...
         'HorizontalAlignment', 'center');
end
ylabel('mV'); title('Envoltório linear e contrações detectadas');
ylim([0 1.15*max(emg_env)]); xlim([0 t(end)]); grid on;

subplot(2,1,2);
plot(t, emg_norm, 'b', 'LineWidth', 1.2);
xlabel('Tempo (s)'); ylabel('% CVIM');
title('RMS normalizado pela CVIM'); xlim([0 t(end)]); grid on;

%% --- Figura: PSD de cada contração com MNF e MDF (FIGURA 8) -------------
figure('Name', 'Domínio da frequência', 'Color', 'w');
cores = lines(nC);
subplot(1,2,1); hold on;
for k = 1:nC
    plot(f, pxx_c(:,k), 'Color', cores(k,:), 'LineWidth', 1.2);
end
xlabel('Frequência (Hz)'); ylabel('PSD (mV^2/Hz)');
title('PSD de cada contração');
legend(compose('C%d', 1:nC)); xlim([0 fs/2]); grid on;

subplot(1,2,2);
k = k_cvim;                       % detalha a contração de referência
plot(f, pxx_c(:,k), 'k', 'LineWidth', 1.2); hold on;
xline(mnf(k), 'b--', sprintf('MNF = %.1f Hz', mnf(k)), 'LineWidth', 1.5);
xline(mdf(k), 'r--', sprintf('MDF = %.1f Hz', mdf(k)), 'LineWidth', 1.5);
xlabel('Frequência (Hz)'); ylabel('PSD (mV^2/Hz)');
title(sprintf('PSD da contração %d (CVIM)', k));
xlim([0 fs/2]); grid on;

%% --- Figura: resumo das medidas por contração ---------------------------
figure('Name', 'Resumo', 'Color', 'w');
subplot(1,3,1);
bar(rms_pico_pct); xlabel('Contração'); ylabel('% CVIM');
title('RMS de pico'); grid on;
subplot(1,3,2);
bar(iemg); xlabel('Contração'); ylabel('mV.s');
title('iEMG'); grid on;
subplot(1,3,3);
bar([mnf mdf]); xlabel('Contração'); ylabel('Hz');
title('Frequências média e mediana'); legend('MNF', 'MDF'); grid on;
