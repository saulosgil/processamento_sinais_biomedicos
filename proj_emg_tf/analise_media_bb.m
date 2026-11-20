close all;
clear all;
clc;

% =========================================================================
%  RMS PICO DA CVIM PARA TODOS OS ARQUIVOS (TEMPO 0 x TEMPO 4)
% =========================================================================

% Canal de EMG a ser analisado
% Escolha o canal desejado: 1 = R BICEPS BRAQUIAL | 21 = R TRICEPS BRAQUIAL
canal   = 1;   

% Pasta com os arquivos .hpf: a mesma pasta deste script (proj_emg_tf)
pasta_dados = fileparts(mfilename('fullpath'));
if isempty(pasta_dados), pasta_dados = pwd; end   % caso o código seja rodado por seção/linha

% Filtros de seleção dos arquivos (sem depender do "ç" de "braço")
braco  = '-t_medida';       % '-t_medida' (treinado) ou '-nt_medida' (não treinado)
medida = 'medida-emg';      % só arquivos de EMG
tarefa = 'tarefa-cvim';
tempos = [0 4];             % coluna 1 = tempo 0, coluna 2 = tempo 4

% Parâmetros do processamento
f_rede   = 60;              % frequência da rede elétrica (Hz)
fc_pb    = [20 450];        % passa-banda (Hz)
ordem_pb = 4;               % ordem do Butterworth passa-banda
t_desc   = 3;               % segundos iniciais descartados (transientes)
janela_s = 0.200;           % janela do RMS móvel (s)

% -------------------------------------------------------------------------
%  1. LISTA E FILTRA OS ARQUIVOS
% -------------------------------------------------------------------------
lista = dir(fullfile(pasta_dados, '*.hpf'));
nomes = {lista.name};
sel   = contains(nomes, tarefa) & contains(nomes, braco) & contains(nomes, medida);
nomes = nomes(sel);

if isempty(nomes)
    error('Nenhum arquivo encontrado em %s com os filtros escolhidos.', pasta_dados);
end

% Extrai sujeito e tempo do nome: "Sub-C001_tempo-0_..."
info     = regexp(nomes, 'Sub-([^_]+)_tempo-(\d+)', 'tokens', 'once');
sujeitos = cellfun(@(x) x{1}, info, 'UniformOutput', false);
tempo_ar = cellfun(@(x) str2double(x{2}), info);

lista_suj = unique(sujeitos);                    % sujeitos em ordem alfabética
nSuj      = numel(lista_suj);

RMS = nan(nSuj, numel(tempos));                  % matriz final: linhas = sujeitos

% -------------------------------------------------------------------------
%  2. LOOP EM TODOS OS ARQUIVOS
% -------------------------------------------------------------------------
for k = 1:numel(nomes)

    col = find(tempos == tempo_ar(k));           % coluna da matriz (tempo 0 ou 4)
    if isempty(col), continue; end               % ignora outros tempos, se houver
    lin = find(strcmp(lista_suj, sujeitos{k}));  % linha da matriz (sujeito)

    arquivo = fullfile(pasta_dados, nomes{k});
    fprintf('[%d/%d] %s\n', k, numel(nomes), nomes{k});

    % ---- Leitura e frequência de amostragem do canal ----
    hpf   = lerHPF(arquivo);
    xml   = strjoin(hpf.xmlCanais, newline);
    taxas = regexp(xml, '<PerChannelSampleRate>([^<]*)</PerChannelSampleRate>', 'tokens');
    taxas = cellfun(@(x) str2double(strrep(x{1}, ',', '.')), taxas);
    fs    = taxas(canal);

    % ---- Pré-processamento ----
    emg_bruto = double(hpf.dados{canal}) * 1e3;  % V -> mV
    emg_bruto = emg_bruto - mean(emg_bruto);     % remove offset DC

    N = numel(emg_bruto);
    t = (0:N-1)' / fs;
    % Remove picos isolados: compara cada amostra com a mediana local
    % (janela de ~50 ms) e substitui as que passam de 6 desvios (MAD)
    jan_out   = round(0.05 * fs);
    emg_bruto = filloutliers(emg_bruto, 'linear', 'movmedian', jan_out, ...
                         'ThresholdFactor', 6);

    % ---- Filtragem ----
    [b_n, a_n]   = butter(2, [f_rede-2 f_rede+2] / (fs/2), 'stop');   % notch 60 Hz
    emg_filt_rede = filtfilt(b_n, a_n, emg_bruto);

    [b_pb, a_pb] = butter(ordem_pb/2, fc_pb / (fs/2), 'bandpass');    % passa-banda de ordem 4
    emg_filt     = filtfilt(b_pb, a_pb, emg_filt_rede);

    % ---- Trecho de análise (exclui os primeiros t_desc segundos) ----
    idx_trecho      = t >= min(t) + t_desc;
    t_trecho        = t(idx_trecho);
    emg_filt_trecho = emg_filt(idx_trecho);

    % ---- RMS móvel e pico ----
    janela  = round(janela_s * fs);
    env_rms = sqrt(movmean(emg_filt_trecho.^2, janela));
    [rms_pico, i_max] = max(env_rms);

    fprintf('      Sujeito %s | tempo %d | RMS pico = %.4f mV em t = %.2f s\n', ...
            sujeitos{k}, tempo_ar(k), rms_pico, t_trecho(i_max));

    % Se houver mais de uma repetição do mesmo sujeito/tempo, guarda a maior
    RMS(lin, col) = max(RMS(lin, col), rms_pico, 'omitnan');
end

% -------------------------------------------------------------------------
%  3. RESULTADO: MATRIZ [tempo 0, tempo 4] E TABELA COM OS SUJEITOS
% -------------------------------------------------------------------------
tabela_rms = array2table(RMS, 'VariableNames', {'Tempo_0', 'Tempo_4'}, ...
                         'RowNames', lista_suj);
disp(tabela_rms)

faltando = any(isnan(RMS), 2);
if any(faltando)
    warning('Sujeitos sem os dois tempos: %s', strjoin(lista_suj(faltando), ', '));
end

% -------------------------------------------------------------------------
%  4. BARPLOT (barras com média e desvio-padrão)
% -------------------------------------------------------------------------
media = mean(RMS, 'omitnan');                 % média de cada coluna (tempo 0 e tempo 4)
dp    = std(RMS, 0, 'omitnan');               % desvio-padrão de cada coluna

figure;
b = bar(1:2, media, 0.8, 'FaceColor', 'flat', 'EdgeColor', 'k');
b.CData = [0.60 0.75 0.90;                    % cor da coluna do tempo 0
           0.95 0.65 0.55];                   % cor da coluna do tempo 4
hold on;

errorbar(1:2, media,zeros(size(dp)), dp, 'k', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 12);   % barras de erro

xticks(1:2); 
xticklabels({'Tempo 0', 'Tempo 4'});
ylabel('RMS pico da CVIM (mV)');
ylim([0 0.24]);
yticks(0:0.05:0.30);
title('RMS pico antes e após o treinamento de força');
grid on;

% -------------------------------------------------------------------------
%  5. Teste-t (pareado) para comparar os tempos 0 e 4
% -------------------------------------------------------------------------
pre = RMS(:, 1);
pos = RMS(:, 2);
dif = pos - pre;

% Teste t para amostras dependentes (pareado): H0 -> média de (pos - pre) = 0
[h, p, ci, stats] = ttest(pos, pre);

% Resultados
fprintf('h = %d  (1 = rejeita H0 a 5%%; 0 = não rejeita)\n', h);
fprintf('p = %.4f\n', p);
fprintf('IC95%% da diferença (pós - pré): %.4f a %.4f mV\n', ci(1), ci(2));
fprintf('t(%d) = %.3f\n', stats.df, stats.tstat);
fprintf('Média das diferenças = %.4f mV\n', mean(dif, 'omitnan'));