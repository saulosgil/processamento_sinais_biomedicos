close all;
clear all;
clc;

% Paramêtros do arquivo de dados
% Canal de EMG a ser analisado.
% Escolha o canal desejado: 1 = R BICEPS BRAQUIAL | 21 = R TRICEPS BRAQUIAL
canal   = 11;           

% Parâmetros do espectro de potência (método de Welch)
janela_psd = 0.5;      % duração de cada janela de Hamming (s)
fmax_psd   = 500;      % frequência máxima exibida no gráfico do espectro (Hz)

%  ------------------------------------------------------------------------
%  1. LISTA DOS ARQUIVOS .HPF DE EMG
%  ------------------------------------------------------------------------
% Procura os arquivos na mesma pasta deste script (independe da pasta atual)
pasta = fileparts(mfilename('fullpath'));
if isempty(pasta), pasta = pwd; end
addpath(pasta);                                  % garante acesso ao lerHPF.m

arquivos = dir(fullfile(pasta, 'Sub-*_medida-emg_*.hpf'));
if isempty(arquivos)
    error('Nenhum arquivo de EMG (.hpf) encontrado em: %s', pasta);
end

% Extrai sujeito (ex.: "Sub-C001") e tempo (ex.: 0, 4) do nome de cada arquivo
nArq    = numel(arquivos);
sujeito = cell(nArq, 1);
tempo   = zeros(nArq, 1);
for k = 1:nArq
    sujeito{k} = regexp(arquivos(k).name, '^Sub-[^_]+', 'match', 'once');
    tempo(k)   = str2double(regexp(arquivos(k).name, 'tempo-(\d+)', 'tokens', 'once'));
end

listaSuj   = unique(sujeito);   % uma linha por sujeito
listaTempo = unique(tempo);     % uma coluna por tempo (0, 4, ...)
nSuj   = numel(listaSuj);
nTempo = numel(listaTempo);

%  ------------------------------------------------------------------------
%  2. LOOP: LEITURA, ESPECTRO DE POTÊNCIA E SINAL BRUTO DE CADA REGISTRO
%  ------------------------------------------------------------------------
% Layout: para cada tempo, 2 colunas -> [Sinal bruto | Espectro]
%   Ex.: Bruto t0 | Espectro t0 | Bruto t4 | Espectro t4
nCol = 2 * nTempo;
figure('Name', 'Sinal EMG bruto e espectro', 'Units', 'normalized', 'Position', [0.02 0.05 0.96 0.85]);
tl = tiledlayout(nSuj, nCol, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:nSuj
    for j = 1:nTempo
        % Posições dos gráficos: linha = sujeito; colunas = bruto e espectro do tempo j
        tile_bruto = (i-1)*nCol + 2*j - 1;
        tile_psd   = (i-1)*nCol + 2*j;
        rotulo = sprintf('%s\\_tempo-%d', listaSuj{i}, listaTempo(j));

        k = find(strcmp(sujeito, listaSuj{i}) & tempo == listaTempo(j), 1);
        if isempty(k)
            % Registro inexistente para esse sujeito/tempo
            for tile = [tile_bruto tile_psd]
                nexttile(tl, tile);
                axis off;
                text(0.5, 0.5, ['Sem registro: ' rotulo], 'HorizontalAlignment', 'center');
            end
            continue;
        end

        arquivo = fullfile(arquivos(k).folder, arquivos(k).name);
        hpf = lerHPF(arquivo);                      % lê todos os canais do arquivo

        % O XML de cada canal traz o nome, a unidade e a frequência de amostragem.
        xml    = strjoin(hpf.xmlCanais, newline);
        nomes  = regexp(xml, '<Name>([^<]*)</Name>', 'tokens');
        nomes  = cellfun(@(x) x{1}, nomes, 'UniformOutput', false);
        taxas  = regexp(xml, '<PerChannelSampleRate>([^<]*)</PerChannelSampleRate>', 'tokens');
        taxas  = cellfun(@(x) str2double(strrep(x{1}, ',', '.')), taxas);
        unids  = regexp(xml, '<Unit>([^<]*)</Unit>', 'tokens');
        unids  = cellfun(@(x) x{1}, unids, 'UniformOutput', false);

        nome_canal = nomes{canal};                  % nome do músculo/canal escolhido
        fs         = taxas(canal);                  % frequência de amostragem (Hz)

        fprintf('%s | Canal %d - %s | Unidade: %s | fs = %.2f Hz\n', ...
                arquivos(k).name, canal, nome_canal, unids{canal}, fs);

        % Sinal de EMG: está em volts, convertido para milivolts (mV)
        emg_bruto = double(hpf.dados{canal}) * 1e3;

        % Remove o nível DC (offset) subtraindo a média do sinal
        emg_bruto = emg_bruto - mean(emg_bruto);
        % Remove picos isolados: compara cada amostra com a mediana local
        
        % (janela de ~50 ms) e substitui as que passam de 6 desvios (MAD)
        jan_out   = round(0.05 * fs);
        emg_bruto = filloutliers(emg_bruto, 'linear', 'movmedian', jan_out, ...
                         'ThresholdFactor', 6);

        % Vetor de tempo (s)
        N = numel(emg_bruto);
        t = (0:N-1)' / fs;

        % Espectro de potência (Welch): janela de Hamming, 50% de sobreposição,
        % nfft = próxima potência de 2 (resolução = fs/nfft)
        nwin  = round(janela_psd * fs);
        nover = round(nwin/2);
        nfft  = 2^nextpow2(nwin);
        [pxx, f] = pwelch(emg_bruto, hamming(nwin), nover, nfft, fs);

        nexttile(tl, tile_psd);
        plot(f, pxx, 'r');
        xlim([0 min(fmax_psd, fs/2)]);
        title(['PSD - ' rotulo]);
        grid on;
        if i == nSuj, xlabel('Frequência (Hz)'); end
        ylabel('PSD (mV^2/Hz)');

        % Sinal Bruto
        nexttile(tl, tile_bruto);
        plot(t, emg_bruto, 'b');
        % ylim([-4 8]);
        % xlim([0 180]);
        title(rotulo);
        grid on;
        if i == nSuj, xlabel('Tempo (s)'); end
        if j == 1,    ylabel('Amplitude (mV)'); end
    end
end

title(tl, ['Sinal EMG bruto e espectro de potência - ' nome_canal]);
