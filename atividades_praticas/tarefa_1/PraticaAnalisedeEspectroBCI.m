% Prática I - Análise Espectral Visando à Construção de Interfaces Cérebro-Máquina

%%%%%%%%%%%%%%%%%%% Exercício 1 %%%%%%%%%%%%%%%%%%%
% Limpa o ambiente e fecha as figuras
clear all;
close all;
clc;

% Carrega os dados
load S03.mat;

%%%%%%%%%%%%%%%%%%% Exercício 2 %%%%%%%%%%%%%%%%%%%
canal = 5;                              % canal POz segundo doc Info.txt
NFFT = 2048;
trials = [1 2];

resolucao = fs/NFFT;                    % resolução espectral (Hz)
fprintf('Resolução espectral: fs/NFFT = %d/%d = %.4f Hz\n', fs, NFFT, resolucao);

figure;
for k = 1:length(trials)
    trial = trials(k);
    sinal = stimuliData(:, canal, trial);
    sinal = sinal - mean(sinal);        % remove o nível DC para não mascarar o espectro
    [Ymag, ~, FreqVect] = EspectAnalysis(sinal, fs, NFFT, 0);

    % Frequência de maior magnitude na faixa dos estímulos (9 a 14 Hz)
    faixa = FreqVect >= 9 & FreqVect <= 14;
    fFaixa = FreqVect(faixa);
    [~, iPico] = max(Ymag(faixa));
    fPico = fFaixa(iPico);

    % Rótulo estimado: estímulo mais próximo do pico
    [~, labelEstimado] = min(abs(frequencies - fPico));
    fprintf('Trial %d: label = %d (%d Hz) | pico em %.3f Hz -> label estimado = %d\n', ...
        trial, label(trial), frequencies(label(trial)), fPico, labelEstimado);

    subplot(length(trials), 1, k);
    plot(FreqVect, Ymag, 'k');
    hold on;
    xline(frequencies, ':b');                                   % frequências de estímulo
    xline(frequencies(label(trial)), 'r', 'LineWidth', 1.5);    % estímulo focado
    hold off;
    xlim([5 30]);
    xlabel('Frequência (Hz)');
    ylabel('|Y(f)|');
    title(sprintf('Canal %d - Trial %d (estímulo de %d Hz)', canal, trial, frequencies(label(trial))));
    grid on;
end

% Respostas:
% - Resolução espectral: fs/NFFT = 256/2048 = 0.125 Hz. Como o sinal tem 1536 amostras (6 s),
%   a resolução "real" é fs/N = 1/6 Hz ~ 0.167 Hz; o zero-padding até 2048 pontos apenas
%   interpola o espectro, sem separar componentes mais próximas que isso.
% - Atribuição do rótulo: nos dois trials o maior pico na faixa 9-14 Hz coincide com o estímulo
%   focado (trial 1: 12 Hz -> label 3; trial 2: 13 Hz -> label 4), então é possível atribuir o rótulo
%   pela frequência de maior magnitude. Porém, no trial 1 o pico em 12 Hz é apenas um pouco maior
%   que os picos vizinhos (~10.2 Hz, ~10.7 Hz) e que a atividade de fundo em baixas frequências: o
%   espectro de um único trial (periodograma simples) tem alta variância, e a decisão é frágil.
%   No trial 2 o pico em 13 Hz se destaca claramente (com harmônico em 26 Hz). Assim, a regra
%   "maior pico -> rótulo" funciona, mas está sujeita a erros; técnicas que reduzem a variância do
%   espectro (ex.: Welch, itens 4-7) tendem a tornar a classificação mais robusta.

%%%%%%%%%%%%%%%%%%% Exercício 3 %%%%%%%%%%%%%%%%%%%
% Função SlidingWindow implementada ao final do arquivo. Teste rápido com um vetor de 10 amostras,
% janelas de 4 pontos e deslocamento de 3 (a última janela é completada com zeros):
[DataWin, inipoint, endpoint, center] = SlidingWindow(1:10, 4, 3)

% Aplicada ao sinal de EEG (1536 amostras, janela de 512, shift de 256): 5 janelas
[DataWin, inipoint, endpoint, center] = SlidingWindow(stimuliData(:, canal, 1), 512, 256);
fprintf('EEG: %d janelas | inipoint = [%s] | endpoint = [%s] | center = [%s]\n', size(DataWin, 2), ...
    num2str(inipoint'), num2str(endpoint'), num2str(center'));

%%%%%%%%%%%%%%%%%%% Exercício 4 %%%%%%%%%%%%%%%%%%%
% Periodograma de Welch (função WelchPeriodogram ao final do arquivo), janela de 512 pontos,
% sobreposição de 50 % (shift = 256) e NFFT = 2048, para o trial 1.
NwindowPoints = 512;
shift = 256;
trial = 1;
sinal = stimuliData(:, canal, trial);
sinal = sinal - mean(sinal);

[PxxHamm, FreqVect] = WelchPeriodogram(sinal, fs, NFFT, NwindowPoints, shift, 'hamming');
[PxxRet, ~]         = WelchPeriodogram(sinal, fs, NFFT, NwindowPoints, shift, 'retangular');
[Ymag, ~, ~]        = EspectAnalysis(sinal, fs, NFFT, 0);

figure;
subplot(3, 1, 1);
plot(FreqVect, Ymag.^2, 'k');
title(sprintf('Trial %d - Espectro de potência do sinal inteiro (|Y|^2, item 1/2)', trial));
subplot(3, 1, 2);
plot(FreqVect, PxxHamm, 'k');
title(sprintf('Trial %d - Welch, janela de Hamming (%d pts, shift %d)', trial, NwindowPoints, shift));
subplot(3, 1, 3);
plot(FreqVect, PxxRet, 'k');
title(sprintf('Trial %d - Welch, janela retangular (%d pts, shift %d)', trial, NwindowPoints, shift));
for k = 1:3
    subplot(3, 1, k);
    hold on;
    xline(frequencies, ':b');
    xline(frequencies(label(trial)), 'r', 'LineWidth', 1.5);
    hold off;
    xlim([5 30]);
    xlabel('Frequência (Hz)');
    ylabel('Potência');
    grid on;
end

% Respostas:
% - Welch x espectro de magnitude do item 1: o periodograma de Welch é bem mais "liso" (menor
%   variância), pois é a média dos espectros de 5 janelas. Os picos espúrios do espectro do sinal
%   inteiro somem e o pico em 12 Hz (estímulo do trial 1) passa a se destacar claramente. O preço é a
%   perda de resolução: cada janela tem 512 pontos (2 s), resolução real de fs/512 = 0.5 Hz, então os
%   picos ficam mais largos.
% - Janela retangular x Hamming: a retangular tem lóbulo principal mais estreito (picos mais finos),
%   porém lóbulos laterais bem mais altos (vazamento espectral, -13 dB contra -43 dB da Hamming); o
%   espectro fica mais "ruidoso", com mais picos pequenos espalhados. Já a Hamming suaviza mais o
%   espectro, mas alarga os picos. Para o trial 1, a retangular parece MELHOR para identificar o
%   estímulo: o pico em 12 Hz (~3.5e-12) fica bem acima dos valores em 10, 11 e 13 Hz (~1.6e-12,
%   ~2.2e-12 e ~0.6e-12). Na Hamming, o lóbulo principal largo funde o pico de ~10.6 Hz com o de
%   12 Hz, e o valor em 11 Hz (~2.6e-12) fica próximo do de 12 Hz (~3.0e-12): o acerto é mais frágil.

%%%%%%%%%%%%%%%%%%% Exercício 5 %%%%%%%%%%%%%%%%%%%
% pwelch do Matlab (usa janela de Hamming por padrão) x implementação do item 4.
[pxx, f] = pwelch(sinal, NwindowPoints, NwindowPoints - shift, NFFT, fs);

figure;
plot(f, pxx/max(pxx), 'k', 'LineWidth', 1.5);
hold on;
plot(FreqVect, PxxHamm/max(PxxHamm), '--r');
xline(frequencies, ':b');
hold off;
xlim([5 30]);
xlabel('Frequência (Hz)');
ylabel('Potência normalizada');
legend('pwelch (Matlab)', 'WelchPeriodogram (item 4)');
title(sprintf('Trial %d - pwelch x implementação própria (Hamming)', trial));
grid on;

R = corrcoef(pxx, PxxHamm);
fprintf('Correlação entre pwelch e o item 4: %.4f\n', R(1, 2));

% Resposta: as duas curvas (normalizadas pelo máximo) se sobrepõem; a diferença é apenas de escala
% (pwelch retorna densidade espectral de potência em V^2/Hz). O pwelch usa o mesmo procedimento:
% janelas de 512 pontos com Hamming, 50 % de sobreposição, média dos periodogramas.

%%%%%%%%%%%%%%%%%%% Exercício 6 %%%%%%%%%%%%%%%%%%%
% Classificação de todos os trials pelo maior valor do periodograma de Welch (Hamming) nas
% frequências de estímulo (função ClassificaSSVEP ao final do arquivo).
[taxaacerto, labelhat] = ClassificaSSVEP(stimuliData(:, canal, :), fs, NFFT, NwindowPoints, ...
    shift, 'hamming', frequencies, label);
fprintf('Canal %d, Hamming, janela %d, NFFT %d: taxa de acerto = %.2f %%\n', ...
    canal, NwindowPoints, NFFT, 100*taxaacerto);
disp(table((1:length(label))', label, labelhat, 'VariableNames', {'trial', 'label', 'labelhat'}));

%%%%%%%%%%%%%%%%%%% Exercício 7 %%%%%%%%%%%%%%%%%%%
% Taxa de acerto para todos os canais, janelas (Hamming e retangular) e NFFT = 256, 512, 1024, 2048.
% A janela tem NFFT pontos e a sobreposição é de 50 % (shift = NFFT/2). Para NFFT = 2048 a janela é
% maior que o sinal (1536 amostras): resta uma única janela completada com zeros.
NFFTs = [256 512 1024 2048];
janelas = {'hamming', 'retangular'};
Ncanais = size(stimuliData, 2);
taxas = zeros(Ncanais, length(janelas), length(NFFTs));     % canal x janela x NFFT

for n = 1:length(NFFTs)
    for j = 1:length(janelas)
        for c = 1:Ncanais
            taxas(c, j, n) = ClassificaSSVEP(stimuliData(:, c, :), fs, NFFTs(n), NFFTs(n), ...
                NFFTs(n)/2, janelas{j}, frequencies, label);
        end
    end
end

% Taxa por canal para a configuração do item 6 (janela de 512 pontos, NFFT = 512)
fprintf('\nTaxa de acerto (%%) por canal, janela = NFFT = 512, shift = 256:\n');
disp(table((1:Ncanais)', 100*taxas(:, 1, 2), 100*taxas(:, 2, 2), ...
    'VariableNames', {'canal', 'hamming', 'retangular'}));

% Desempenho de pico (melhor canal) para cada janela e NFFT
fprintf('Desempenho de pico (melhor canal) por janela e NFFT:\n');
for n = 1:length(NFFTs)
    for j = 1:length(janelas)
        [pico, melhorCanal] = max(taxas(:, j, n));
        Nwin = size(SlidingWindow(sinal, NFFTs(n), NFFTs(n)/2), 2);
        fprintf('NFFT = %4d (%d janelas, df = %.3f Hz) | %-10s | pico = %6.2f %% (canal %d)\n', ...
            NFFTs(n), Nwin, fs/NFFTs(n), janelas{j}, 100*pico, melhorCanal);
    end
end

figure;
plot(NFFTs, 100*squeeze(max(taxas(:, 1, :), [], 1)), '-ok', 'LineWidth', 1.5);
hold on;
plot(NFFTs, 100*squeeze(max(taxas(:, 2, :), [], 1)), '--sr', 'LineWidth', 1.5);
hold off;
set(gca, 'XScale', 'log', 'XTick', NFFTs);
xlabel('NFFT (= tamanho da janela)');
ylabel('Taxa de acerto de pico (%)');
legend('Hamming', 'Retangular', 'Location', 'best');
title('Desempenho de pico entre os canais');
grid on;

figure;
imagesc(100*squeeze(max(taxas, [], 2)));
colorbar;
set(gca, 'XTick', 1:length(NFFTs), 'XTickLabel', NFFTs, 'YTick', 1:Ncanais);
xlabel('NFFT');
ylabel('Canal');
title('Taxa de acerto (%) por canal (melhor janela)');

% Respostas (valores obtidos na execução):
% a) Melhor canal: o canal 1 leva ao melhor desempenho em todas as configurações (pico de 72.92 %),
%    seguido dos canais 3 e 4. (Associar o índice ao nome do eletrodo requer o Info.txt.)
% b) Melhor janela: a retangular, com desempenho de pico maior ou igual ao da Hamming em todos os
%    NFFT (72.92 % x 70.83 % no melhor caso). Com poucas janelas, a resolução (lóbulo principal
%    estreito) pesa mais que o vazamento, pois os estímulos estão separados por apenas 1 Hz.
% c) Desempenho de pico (Hamming / retangular):
%      NFFT =  256 -> 11 janelas, df = 1.000 Hz: 66.67 % / 70.83 %
%      NFFT =  512 ->  5 janelas, df = 0.500 Hz: 66.67 % / 72.92 %
%      NFFT = 1024 ->  2 janelas, df = 0.250 Hz: 70.83 % / 72.92 %
%      NFFT = 2048 ->  1 janela,  df = 0.125 Hz: 70.83 % / 72.92 %
%    Compromisso número de janelas x resolução: janelas curtas permitem mais janelas na média (menor
%    variância do espectro), mas pioram a resolução. Com 256 pontos, df = 1 Hz é igual ao espaçamento
%    entre os estímulos, e o lóbulo principal de um estímulo invade as frequências vizinhas. Janelas
%    longas melhoram a resolução, mas sobram poucas janelas para a média (com 2048 pontos há uma única
%    janela com zero-padding, i.e., o periodograma simples, cuja resolução real é limitada pelos 6 s do
%    sinal: 1/6 Hz). Aqui a resolução pesa mais: o desempenho cresce de 256 para 1024 e satura a partir
%    daí, já que acima disso o ganho de resolução é só interpolação e não há mais janelas para reduzir a
%    variância.

% Exercício 1
function [Ymag, Yangle, FreqVect] = EspectAnalysis(Sinal, Fs, NFFT, grafico)
% EspectAnalysis calcula o espectro de magnitude de lado único (single-sided).
%   [Ymag, Yangle, FreqVect] = EspectAnalysis(Sinal, Fs, NFFT, grafico)
%   Sinal   - sinal no tempo (vetor)
%   Fs      - frequência de amostragem (Hz)
%   NFFT    - número de pontos da FFT
%   grafico - 1 (true) para plotar o espectro de magnitude, 0 (false) caso contrário
%   Ymag     - magnitude do espectro de lado único, de 0 a Fs/2
%   Yangle   - fase (rad) correspondente
%   FreqVect - vetor de frequências (Hz)

if nargin < 4
    grafico = 0;
end

Sinal = Sinal(:);                       % garante vetor coluna
N = length(Sinal);                      % número de amostras do sinal

Y = fft(Sinal, NFFT)/N;                 % double-sided spectrum normalizado pelo número de amostras
Nhalf = floor(NFFT/2) + 1;              % coeficientes de 0 a Fs/2

Ymag = abs(Y(1:Nhalf));
Ymag(2:end) = 2*Ymag(2:end);            % dobra a magnitude em (0, Fs/2], exceto o DC (f = 0)
Yangle = angle(Y(1:Nhalf));
FreqVect = (0:Nhalf-1)'*Fs/NFFT;        % resolução espectral = Fs/NFFT

if grafico
    figure;
    plot(FreqVect, Ymag, 'k');
    xlabel('Frequência (Hz)');
    ylabel('|Y(f)|');
    title('Espectro de magnitude (single-sided)');
    grid on;
end
end

% Exercício 3
function [DataWin, inipoint, endpoint, center] = SlidingWindow(Data, NwindowPoints, shift)
% SlidingWindow fatia o sinal em janelas deslizantes.
%   [DataWin, inipoint, endpoint, center] = SlidingWindow(Data, NwindowPoints, shift)
%   Data          - vetor no tempo
%   NwindowPoints - número de amostras de cada janela
%   shift         - deslocamento entre o início de janelas consecutivas
%   DataWin  - matriz NwindowPoints x N janelas (a última é completada com zeros, se necessário)
%   inipoint - índice de início de cada janela
%   endpoint - índice de fim de cada janela
%   center   - índice do centro de cada janela, ceil(median(inipoint:endpoint))

Data = Data(:);
L = length(Data);

% Índices de início: a última janela é a primeira que alcança o fim do sinal
inipoint = 1;
while inipoint(end) + NwindowPoints - 1 < L
    inipoint(end+1) = inipoint(end) + shift;
end
inipoint = inipoint(:);
endpoint = inipoint + NwindowPoints - 1;

Nwin = length(inipoint);
DataWin = zeros(NwindowPoints, Nwin);
center = zeros(Nwin, 1);
for k = 1:Nwin
    idx = inipoint(k):min(endpoint(k), L);              % índices disponíveis no sinal
    DataWin(1:length(idx), k) = Data(idx);              % restante fica com zeros
    center(k) = ceil(median(inipoint(k):endpoint(k)));
end
end

% Exercício 4
function [Pxx, FreqVect] = WelchPeriodogram(Sinal, Fs, NFFT, NwindowPoints, shift, tipoJanela)
% WelchPeriodogram calcula o espectro de potência pelo periodograma de Welch.
%   [Pxx, FreqVect] = WelchPeriodogram(Sinal, Fs, NFFT, NwindowPoints, shift, tipoJanela)
%   tipoJanela - 'hamming' ou 'retangular'

DataWin = SlidingWindow(Sinal, NwindowPoints, shift);

if strcmpi(tipoJanela, 'hamming')
    w = hamming(NwindowPoints);
else
    w = ones(NwindowPoints, 1);
end
U = mean(w.^2);                         % compensa a perda de potência causada pela janela

Pxx = 0;
for k = 1:size(DataWin, 2)
    [Ymag, ~, FreqVect] = EspectAnalysis(DataWin(:, k).*w, Fs, NFFT, 0);
    Pxx = Pxx + Ymag.^2;
end
Pxx = Pxx/(size(DataWin, 2)*U);         % média dos periodogramas das janelas
end

% Exercícios 6 e 7
function [taxaacerto, labelhat] = ClassificaSSVEP(Dados, Fs, NFFT, NwindowPoints, shift, tipoJanela, frequencies, label)
% ClassificaSSVEP atribui a cada trial a classe cuja frequência de estímulo tem a maior potência no
% periodograma de Welch (máximo em +-0.25 Hz ao redor de cada frequência) e calcula a taxa de acerto.
%   Dados - amostras x 1 canal x trials

Dados = squeeze(Dados);
Ntrials = size(Dados, 2);
labelhat = zeros(Ntrials, 1);
for k = 1:Ntrials
    sinal = Dados(:, k) - mean(Dados(:, k));
    [Pxx, FreqVect] = WelchPeriodogram(sinal, Fs, NFFT, NwindowPoints, shift, tipoJanela);
    potencia = zeros(1, length(frequencies));
    for i = 1:length(frequencies)
        potencia(i) = max(Pxx(abs(FreqVect - frequencies(i)) <= 0.25));
    end
    [~, labelhat(k)] = max(potencia);
end
taxaacerto = 1 - (sum(labelhat ~= label)/length(label));
end
