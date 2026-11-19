
%% Análise comparativa de EMG: saudável x miopatia x neuropatia
% Dados: PhysioNet EMGDB (https://physionet.org/content/emgdb/1.0.0/)
% EMG de agulha do tibial anterior, fs = 4 kHz, amplitude em mV.
%
% Etapas:
%   1. Leitura robusta dos dados
%   2. Pré-processamento (verificação de rede elétrica + passa-banda IIR x FIR)
%   3. Detecção e média dos potenciais de unidade motora (MUAPs)
%   4. Características por janela (tempo, frequência e análise turns/amplitude)
%   5. Espectrograma (STFT) e escalograma (wavelet de Morlet)
%   6. Separação dos grupos (dispersão, PCA e classificador kNN)
%
% Requer apenas: MATLAB + Signal Processing Toolbox.

clear; close all; clc;

%% Configurações gerais
pasta    = fileparts(mfilename('fullpath'));
arquivos = {'emg_healthy.txt', 'emg_myopathy.txt', 'emg_neuropathy.txt'};
titulos  = {'Saudável', 'Miopatia', 'Neuropatia'};
cores    = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880]};
nS       = numel(arquivos);

fs            = 4000;      % Hz (informado na documentação do EMGDB)
banda         = [20 1500]; % Hz, passa-banda
tJanela       = 0.5;       % s, janela para extração de características
limiarTurn    = 0.1;       % mV (100 µV, critério clássico de Willison)
salvarFiguras = true;      % salva PNGs na pasta "figuras" (para os slides)

pastaFig = fullfile(pasta, 'figuras');
if salvarFiguras && ~exist(pastaFig, 'dir'), mkdir(pastaFig); end

%% 1. Leitura dos dados
% Os arquivos têm algumas linhas corrompidas (caracteres inválidos) e saltos
% na coluna de tempo. Por isso: lê-se apenas o 1º valor numérico da coluna de
% amplitude e o eixo de tempo é reconstruído a partir de fs.
x = cell(1, nS);  t = cell(1, nS);
fprintf('--- Leitura ---\n');
for k = 1:nS
    [x{k}, nCorr] = ler_emg(fullfile(pasta, arquivos{k}));
    t{k} = (0:numel(x{k})-1)' / fs;
    fprintf('%-10s | N = %6d | duração = %5.2f s | linhas corrigidas = %d\n', ...
        titulos{k}, numel(x{k}), t{k}(end), nCorr);
end

%% 2. Pré-processamento
% 2a. Verifica interferência da rede elétrica (50/60 Hz) na PSD
fprintf('\n--- Interferência da rede (pico acima da vizinhança) ---\n');
for k = 1:nS
    [P, f] = pwelch(x{k} - mean(x{k}), hamming(fs), fs/2, fs, fs);  % resolução 1 Hz
    for fr = [50 60]
        viz = abs(f - fr) <= 10 & abs(f - fr) > 2;
        ganho = 10*log10(P(f == fr) / median(P(viz)));
        fprintf('%-10s | %d Hz: %+5.1f dB\n', titulos{k}, fr, ganho);
        if ganho > 6   % notch apenas se houver pico evidente
            [bn, an] = iirnotch(fr/(fs/2), (fr/(fs/2))/35);
            x{k} = filtfilt(bn, an, x{k});
            fprintf('           -> notch aplicado em %d Hz\n', fr);
        end
    end
end

% 2b. Projeto dos filtros passa-banda: IIR (Butterworth) x FIR (janela de Hamming)
[z, p, g]  = butter(4, banda/(fs/2), 'bandpass');    % ordem 8
[sos, gS]  = zp2sos(z, p, g);
bFIR       = fir1(400, banda/(fs/2), 'bandpass');    % ordem 400

nF = 4096;
[hI, fI] = freqz(sos, nF, fs);   hI = hI * gS;
[hF, fF] = freqz(bFIR, 1, nF, fs);

fig = figure('Name', 'Filtros', 'Color', 'w');
subplot(2,1,1);
plot(fI, 20*log10(abs(hI)), fF, 20*log10(abs(hF)), 'LineWidth', 1.2);
grid on; ylim([-80 5]); xlim([0 fs/2]);
title('Resposta em magnitude'); ylabel('|H(f)| (dB)');
legend('IIR Butterworth (ordem 8)', 'FIR Hamming (ordem 400)', 'Location', 'south');
subplot(2,1,2);
plot(fI, unwrap(angle(hI)), fF, unwrap(angle(hF)), 'LineWidth', 1.2);
grid on; xlim([0 fs/2]);
title('Resposta em fase (IIR: não linear | FIR: linear)');
xlabel('Frequência (Hz)'); ylabel('Fase (rad)');
salvar(fig, pastaFig, '01_filtros', salvarFiguras);

% Aplica o IIR com filtfilt (fase zero: não distorce a forma dos MUAPs)
xf = cell(1, nS);
for k = 1:nS
    xf{k} = filtfilt(sos, gS, x{k} - mean(x{k}));
end

% Efeito do filtro: filter (fase não linear) x filtfilt (fase zero)
idx = t{1} >= 1 & t{1} <= 1.1;
xFilter = gS * sosfilt(sos, x{1} - mean(x{1}));   % filtragem causal (uma passada)
fig = figure('Name', 'Efeito do filtro', 'Color', 'w');
plot(t{1}(idx), x{1}(idx) - mean(x{1}), 'Color', [0.6 0.6 0.6]); hold on;
plot(t{1}(idx), xFilter(idx), 'r');
plot(t{1}(idx), xf{1}(idx), 'b', 'LineWidth', 1.2);
grid on; xlabel('Tempo (s)'); ylabel('Amplitude (mV)');
title('Sinal saudável: bruto x filtrado');
legend('Bruto', 'sosfilt - causal (distorção de fase)', 'filtfilt (fase zero)');
salvar(fig, pastaFig, '02_efeito_filtro', salvarFiguras);

%% 3. Detecção e média dos MUAPs
% Detecta picos acima de 5 desvios-padrão do ruído (estimativa robusta via
% MAD), recorta ±15 ms em torno de cada pico, alinha pelo pico e corrige a
% polaridade (pico principal positivo).
meia   = round(0.015 * fs);           % 15 ms
tMUAP  = (-meia:meia)' / fs * 1000;   % ms
MU     = struct('ondas', [], 'Vpp', [], 'dur', [], 'fases', []);
fprintf('\n--- MUAPs ---\n');
for k = 1:nS
    s     = xf{k};
    sigma = median(abs(s)) / 0.6745;
    [~, loc] = findpeaks(abs(s), 'MinPeakHeight', 5*sigma, ...
                         'MinPeakDistance', round(0.010*fs));
    loc = loc(loc > meia & loc <= numel(s) - meia);

    ondas = zeros(numel(tMUAP), numel(loc));
    for i = 1:numel(loc)
        seg = s(loc(i)-meia : loc(i)+meia);
        ondas(:, i) = seg * sign(seg(meia+1));
    end

    [Vpp, dur, fases] = medidas_muap(ondas, fs, sigma);
    MU(k).ondas = ondas;  MU(k).Vpp = Vpp;  MU(k).dur = dur;  MU(k).fases = fases;

    fprintf(['%-10s | n = %4d | Vpp = %.2f mV | duração = %.1f ms | ' ...
             'fases = %.1f | polifásicos (>4 fases) = %.0f%%\n'], ...
        titulos{k}, numel(loc), median(Vpp), median(dur), median(fases), ...
        100*mean(fases > 4));
end

% Figura: ondas sobrepostas + mediana
fig = figure('Name', 'MUAPs', 'Color', 'w');
yl = max(cellfun(@(o) max(abs(median(o, 2))), {MU.ondas})) * 1.3;  % mesma escala p/ comparar
for k = 1:nS
    subplot(1, nS, k);
    o = MU(k).ondas;
    sel = round(linspace(1, size(o, 2), min(150, size(o, 2))));
    plot(tMUAP, o(:, sel), 'Color', [0.8 0.8 0.8]); hold on;
    plot(tMUAP, median(o, 2), 'Color', cores{k}, 'LineWidth', 2.5);
    grid on; ylim([-yl yl]);
    title(sprintf('%s (n = %d)', titulos{k}, size(o, 2)));
    xlabel('Tempo (ms)');
    if k == 1, ylabel('Amplitude (mV)'); end
end
sgtitle('MUAPs detectados (cinza) e forma mediana (cor)');
salvar(fig, pastaFig, '03_muaps', salvarFiguras);

% Figura: distribuição das medidas dos MUAPs
fig = figure('Name', 'Medidas dos MUAPs', 'Color', 'w');
med = {'Vpp', 'dur', 'fases'};
rot = {'Amplitude pico a pico (mV)', 'Duração (ms)', 'Número de fases'};
for m = 1:3
    subplot(1, 3, m);
    [v, gcat] = empilha(cellfun(@(c) c(:), {MU.(med{m})}, 'UniformOutput', false), titulos);
    boxchart(gcat, v, 'MarkerStyle', 'none'); grid on;
    title(rot{m});
end
salvar(fig, pastaFig, '04_medidas_muaps', salvarFiguras);

%% 4. Características por janela
nomesF = {'RMS (mV)', 'MAV (mV)', 'Vpp (mV)', 'Cruz. zero (/s)', ...
          'Comprimento de onda (mV/s)', 'Curtose', ...
          'f média (Hz)', 'f mediana (Hz)', ...
          'P 20-150 Hz (%)', 'P 150-500 Hz (%)', 'P 500-1500 Hz (%)', ...
          'Turns (/s)', 'Amplitude média/turn (mV)'};
nJ = round(tJanela * fs);
F = [];  grupo = [];  ordem = [];
for k = 1:nS
    nW = floor(numel(xf{k}) / nJ);
    for w = 1:nW
        seg = xf{k}((w-1)*nJ + 1 : w*nJ);
        F(end+1, :) = caracteristicas(seg, fs, limiarTurn); %#ok<SAGROW>
        grupo(end+1, 1) = k;                                %#ok<SAGROW>
        ordem(end+1, 1) = w / nW;                           %#ok<SAGROW>  posição relativa no tempo
    end
end
gCat = categorical(grupo, 1:nS, titulos);

fprintf('\n--- Características (mediana das janelas de %.1f s) ---\n', tJanela);
fprintf('%-28s', ''); fprintf('%12s', titulos{:}); fprintf('\n');
for c = 1:numel(nomesF)
    fprintf('%-28s', nomesF{c});
    for k = 1:nS, fprintf('%12.2f', median(F(grupo == k, c))); end
    fprintf('\n');
end

fig = figure('Name', 'Características por janela', 'Color', 'w', ...
             'Position', [50 50 1300 800]);
tl = tiledlayout(4, 4, 'TileSpacing', 'compact');
for c = 1:numel(nomesF)
    nexttile;
    boxchart(gCat, F(:, c), 'MarkerStyle', '.'); grid on;
    title(nomesF{c}, 'FontWeight', 'normal');
end
title(tl, sprintf('Características em janelas de %.1f s', tJanela));
salvar(fig, pastaFig, '05_caracteristicas', salvarFiguras);

% Potência relativa por banda (barras empilhadas)
fig = figure('Name', 'Potência por banda', 'Color', 'w');
Pb = zeros(nS, 3);
for k = 1:nS, Pb(k, :) = median(F(grupo == k, 9:11), 1); end
Pb = 100 * Pb ./ sum(Pb, 2);
bar(categorical(titulos, titulos), Pb, 'stacked'); grid on;
ylabel('Potência relativa (%)');
legend('20-150 Hz', '150-500 Hz', '500-1500 Hz', 'Location', 'eastoutside');
title('Distribuição da potência por banda de frequência');
salvar(fig, pastaFig, '06_potencia_bandas', salvarFiguras);

%% 5. Tempo-frequência
% 5a. Espectrograma (STFT) de 3 s de cada sinal
tSeg = [1 4];   % s
fig = figure('Name', 'Espectrogramas', 'Color', 'w');
for k = 1:nS
    subplot(nS, 1, k);
    seg = xf{k}(t{k} >= tSeg(1) & t{k} < tSeg(2));
    [~, fS, tS, pS] = spectrogram(seg, hamming(256), 192, 512, fs);
    pdB = 10*log10(pS);
    imagesc(tS + tSeg(1), fS, pdB); axis xy;
    clim(max(pdB(:)) + [-50 0]); colorbar;
    ylim(banda); ylabel('f (Hz)');
    title(['Espectrograma - ' titulos{k}]);
end
xlabel('Tempo (s)'); colormap(turbo);
salvar(fig, pastaFig, '07_espectrogramas', salvarFiguras);

% 5b. Escalograma com wavelet de Morlet (implementação própria) em 0,5 s
tSeg = [1 1.5];
fW   = logspace(log10(banda(1)), log10(banda(2)), 60);
fig  = figure('Name', 'Escalogramas', 'Color', 'w');
for k = 1:nS
    i0  = t{k} >= tSeg(1) & t{k} < tSeg(2);
    seg = xf{k}(i0);  tt = t{k}(i0);
    W   = cwt_morlet(seg, fs, fW);

    subplot(nS, 2, 2*k - 1);
    plot(tt, seg, 'Color', cores{k}); grid on; xlim(tSeg); ylabel('mV');
    title(['Sinal - ' titulos{k}]);
    subplot(nS, 2, 2*k);
    surface(tt, fW, abs(W), 'EdgeColor', 'none');
    set(gca, 'YScale', 'log'); axis tight; view(2); ylabel('f (Hz)');
    title(['Escalograma - ' titulos{k}]);
end
subplot(nS, 2, 2*nS - 1); xlabel('Tempo (s)');
subplot(nS, 2, 2*nS);     xlabel('Tempo (s)'); colormap(turbo);
salvar(fig, pastaFig, '08_escalogramas', salvarFiguras);

%% 6. Separação dos grupos
% 6a. Dispersão de pares de características
pares = {[1 8], [12 13]};   % RMS x f mediana | turns x amplitude/turn
fig = figure('Name', 'Dispersão', 'Color', 'w', 'Position', [100 100 1100 450]);
for p = 1:numel(pares)
    subplot(1, numel(pares), p); hold on;
    for k = 1:nS
        sel = grupo == k;
        scatter(F(sel, pares{p}(1)), F(sel, pares{p}(2)), 30, cores{k}, 'filled', ...
            'MarkerFaceAlpha', 0.7);
    end
    grid on; xlabel(nomesF{pares{p}(1)}); ylabel(nomesF{pares{p}(2)});
    legend(titulos, 'Location', 'best');
end
subplot(1, 2, 1); set(gca, 'XScale', 'log'); title('Amplitude x conteúdo espectral');
subplot(1, 2, 2); title('Análise turns/amplitude');
salvar(fig, pastaFig, '09_dispersao', salvarFiguras);

% 6b. PCA (via SVD das características padronizadas)
Z = (F - mean(F)) ./ std(F);
[~, S, V] = svd(Z, 'econ');
scores    = Z * V;
explic    = 100 * diag(S).^2 / sum(diag(S).^2);
fig = figure('Name', 'PCA', 'Color', 'w');
hold on;
for k = 1:nS
    sel = grupo == k;
    scatter(scores(sel, 1), scores(sel, 2), 30, cores{k}, 'filled', 'MarkerFaceAlpha', 0.7);
end
grid on; legend(titulos, 'Location', 'best');
xlabel(sprintf('CP1 (%.0f%% da variância)', explic(1)));
ylabel(sprintf('CP2 (%.0f%% da variância)', explic(2)));
title('PCA das características por janela');
salvar(fig, pastaFig, '10_pca', salvarFiguras);

% 6c. Classificador kNN (k = 5) com divisão TEMPORAL:
% treino = 1ª metade de cada sinal, teste = 2ª metade.
% (evita que janelas vizinhas, quase idênticas, caiam no treino e no teste)
treino = ordem <= 0.5;  teste = ~treino;
mu = mean(F(treino, :));  sd = std(F(treino, :));
Ztr = (F(treino, :) - mu) ./ sd;   ytr = grupo(treino);
Zte = (F(teste, :)  - mu) ./ sd;   yte = grupo(teste);
D   = sum(Zte.^2, 2) + sum(Ztr.^2, 2)' - 2 * (Zte * Ztr');
[~, viz] = sort(D, 2);
yPred = mode(ytr(viz(:, 1:5)), 2);

C = accumarray([yte yPred], 1, [nS nS]);
fprintf('\n--- kNN (k = 5), divisão temporal ---\n');
fprintf('Acurácia no teste: %.1f%% (%d janelas)\n', 100*mean(yPred == yte), numel(yte));
fig = figure('Name', 'Matriz de confusão', 'Color', 'w');
heatmap(titulos, titulos, C, 'Colormap', parula);
xlabel('Classe prevista'); ylabel('Classe real');
title(sprintf('kNN (k = 5) - acurácia = %.1f%%', 100*mean(yPred == yte)));
salvar(fig, pastaFig, '11_confusao', salvarFiguras);

fprintf(['\nATENÇÃO: há apenas 1 indivíduo por grupo. O classificador pode estar\n' ...
         'separando INDIVÍDUOS, não DOENÇAS. Use os resultados como exploração\n' ...
         'das características, não como validação diagnóstica.\n']);


%% ======================== Funções locais ========================

function [x, nCorr] = ler_emg(arquivo)
% Lê o 2º valor numérico de cada linha (amplitude), tolerando lixo no final.
    txt    = fileread(arquivo);
    linhas = regexp(txt, '\r?\n', 'split');
    linhas = linhas(~cellfun(@isempty, strtrim(linhas)));
    num    = '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?';
    tok    = regexp(linhas, ['^\s*' num '\s+(' num ')'], 'tokens', 'once');
    ok     = ~cellfun(@isempty, tok);
    tok    = [tok{ok}];
    x      = str2double(tok(:));
    limpas = ~cellfun(@isempty, regexp(linhas, ['^\s*' num '\s+' num '\s*$'], 'once'));
    nCorr  = sum(~limpas);
end

function [Vpp, dur, fases] = medidas_muap(ondas, fs, sigma)
% Amplitude pico a pico, duração e número de fases de cada MUAP.
% Duração: região contínua em torno do pico em que a envoltória (|x| suavizado
% em 1 ms) fica acima de max(10% do pico, 2*sigma do ruído).
    n = size(ondas, 2);  c = (size(ondas, 1) + 1) / 2;
    Vpp = zeros(1, n);  dur = zeros(1, n);  fases = zeros(1, n);
    for i = 1:n
        o   = ondas(:, i);
        Vpp(i) = max(o) - min(o);
        env = movmean(abs(o), round(0.001*fs));
        acima = env > max(0.1*max(abs(o)), 2*sigma);
        ini = find(~acima(1:c), 1, 'last');  if isempty(ini), ini = 0; end
        fim = find(~acima(c:end), 1, 'first'); if isempty(fim), fim = numel(o) - c + 2; end
        ini = ini + 1;  fim = c + fim - 2;
        dur(i) = (fim - ini + 1) / fs * 1000;
        r = o(ini:fim);
        r = r(abs(r) > sigma);                 % ignora oscilações dentro do ruído
        fases(i) = sum(diff(sign(r)) ~= 0) + 1;
    end
end

function v = caracteristicas(s, fs, limiarTurn)
% Vetor de características de uma janela de sinal.
    N    = numel(s);  T = N / fs;
    rmsv = sqrt(mean(s.^2));
    mav  = mean(abs(s));
    vpp  = max(s) - min(s);
    zc   = sum(s(1:end-1) .* s(2:end) < 0 & abs(diff(s)) >= 0.02) / T;  % zona morta 20 µV
    wl   = sum(abs(diff(s))) / T;
    m    = s - mean(s);
    kurt = mean(m.^4) / mean(m.^2)^2;

    [P, f] = pwelch(s, hamming(512), 256, 512, fs);
    b      = f >= 20 & f <= 1500;  P = P(b);  f = f(b);
    fmed   = sum(f .* P) / sum(P);
    Pc     = cumsum(P) / sum(P);
    fmd    = f(find(Pc >= 0.5, 1));
    pb     = 100 * [sum(P(f < 150)), sum(P(f >= 150 & f < 500)), sum(P(f >= 500))] / sum(P);

    [nT, ampT] = turns(s, limiarTurn);
    v = [rmsv, mav, vpp, zc, wl, kurt, fmed, fmd, pb, nT / T, ampT];
end

function [n, ampMed] = turns(s, limiar)
% Análise turns/amplitude (Willison): conta inversões de direção do sinal
% com variação de amplitude >= limiar em relação ao turn anterior.
    ext = find(diff(sign(diff(s))) ~= 0) + 1;
    if isempty(ext), n = 0; ampMed = 0; return; end
    ult = s(ext(1));  amps = zeros(numel(ext), 1);  n = 0;
    for i = 2:numel(ext)
        d = s(ext(i)) - ult;
        if abs(d) >= limiar
            n = n + 1;  amps(n) = abs(d);  ult = s(ext(i));
        end
    end
    if n > 0, ampMed = mean(amps(1:n)); else, ampMed = 0; end
end

function W = cwt_morlet(s, fs, freqs)
% Transformada wavelet contínua com wavelet de Morlet complexa (w0 = 6).
    w0 = 6;  s = s(:);  W = zeros(numel(freqs), numel(s));
    for i = 1:numel(freqs)
        sig = w0 / (2*pi*freqs(i));                  % desvio no tempo (s)
        tw  = (-ceil(4*sig*fs):ceil(4*sig*fs))' / fs;
        psi = exp(1i*2*pi*freqs(i)*tw) .* exp(-tw.^2 / (2*sig^2));
        psi = psi / sum(abs(psi));
        W(i, :) = conv(s, conj(flipud(psi)), 'same');
    end
end

function [v, g] = empilha(c, nomes)
% Junta células de vetores em um vetor e um rótulo categórico.
    v = vertcat(c{:});
    g = categorical(repelem(1:numel(c), cellfun(@numel, c))', 1:numel(c), nomes);
end

function salvar(fig, pasta, nome, flag)
    if flag, exportgraphics(fig, fullfile(pasta, [nome '.png']), 'Resolution', 200); end
end
