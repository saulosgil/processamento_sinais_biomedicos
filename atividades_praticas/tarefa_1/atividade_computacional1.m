% Prática I - Análise Espectral Visando à Construção de Interfaces Cérebro-Máquina
% Saulo Gil
% 24/10/2026
%
close all;
clear all;
clc;

% Carrega os dados
load S03.mat;

% Parâmetros gerais
canal = 5;          % canal POz (segundo o arquivo Info.txt)
NFFT = 2048;        % número de pontos da FFT

%%%%%%%%%%%%%%%%%%% Exercício 1 %%%%%%%%%%%%%%%%%%%
% A função EspectAnalysis 
function [Ymag, Yangle, FreqVect] = EspectAnalysis(Sinal, Fs, NFFT, Plotagrafico)
% Argumentos da função EspectAnalysis (enunciado do exercício).
%   Sinal         - sinal no tempo (vetor)
%   Fs            - frequência de amostragem (Hz)
%   NFFT          - número de pontos da FFT
%   Plotagrafico  - 1 (true) para plotar, 0 (false) para não plotar
% Resultados
%   Ymag          - magnitude do espectro de lado único, de 0 a Fs/2
%   Yangle        - fase (rad) correspondente
%   FreqVect      - vetor de frequências (Hz)

if nargin < 4
    Plotagrafico = 0;
end

N = length(Sinal);                      % número de amostras do sinal

Y = fft(Sinal, NFFT)/N;                 % double-sided spectrum normalizado pelo número de amostras
Nhalf = floor(NFFT/2) + 1;              % coeficientes de 0 a Fs/2

Ymag = abs(Y(1:Nhalf));
Ymag = Ymag(:);
Ymag(2:end) = 2*Ymag(2:end);            % dobra a magnitude em (0, Fs/2], exceto o DC (f = 0)
Yangle = angle(Y(1:Nhalf));
Yangle = Yangle(:);
FreqVect = (0:Nhalf-1)'*Fs/NFFT;        % resolução espectral = Fs/NFFT

if Plotagrafico
    figure;
    plot(FreqVect, Ymag, 'k');
    xlabel('Frequência (Hz)');
    ylabel('|Y(f)|');
    title('Espectro de magnitude');
    grid on;
end
end

%%%%%%%%%%%%%%%%%%% Exercício 2 %%%%%%%%%%%%%%%%%%%
% Resolução do espectro
resolucao = fs/NFFT;
disp(['Qual a resolução espectral?']);
disp(['Resolução espectral = ' num2str(resolucao) ' Hz']);
disp(['-----------------------------------------------------']);

% Realização 1
sinal1 = stimuliData(:, canal, 1);
[Ymag1, fase1, f1] = EspectAnalysis(sinal1, fs, NFFT, 0);

% Pico do espectro entre 9 e 14 Hz (faixa dos estímulos)
faixa1 = find(f1 >= 9 & f1 <= 14);
[valor1, posicao1] = max(Ymag1(faixa1));
pico1 = f1(faixa1(posicao1));

disp(['Realização 1: estímulo focado = ' num2str(frequencies(label(1))) ' Hz']);
disp(['Realização 1: pico do espectro = ' num2str(pico1) ' Hz']);
disp(['-----------------------------------------------------']);

% Realização 2
sinal2 = stimuliData(:, canal, 2);
[Ymag2, fase2, f2] = EspectAnalysis(sinal2, fs, NFFT, 0);

% Pico do espectro entre 9 e 14 Hz (faixa dos estímulos)
faixa2 = find(f2 >= 9 & f2 <= 14);
[valor2, posicao2] = max(Ymag2(faixa2));
pico2 = f2(faixa2(posicao2));

disp(['Realização 2: estímulo focado = ' num2str(frequencies(label(2))) ' Hz']);
disp(['Realização 2: pico do espectro = ' num2str(pico2) ' Hz']);
disp(['-----------------------------------------------------']);

% Gráficos
figure;

subplot(2,1,1);
plot(f1, Ymag1, 'k');
xlim([5 30]);
xlabel('Frequência (Hz)');
ylabel('|Y(f)|');
title(['POz - Realização 1 - Estímulo de ' num2str(frequencies(label(1))) ' Hz']);
grid on;

subplot(2,1,2);
plot(f2, Ymag2, 'k');
xlim([5 30]);
xlabel('Frequência (Hz)');
ylabel('|Y(f)|');
title(['POz - Realização 2 - Estímulo de ' num2str(frequencies(label(2))) ' Hz']);
grid on;

% Resposta Final da 2
% A análise dos resultados indicam que que atribuir o rótulo apenas pelo maior pico do espectro 
% nem sempre funciona. Das duas realizações analisadas, a tentativa 1 não foi classificada corretamente,
% já a tentativa 2 foi classificada corretamente.

%%%%%%%%%%%%%%%%%%% Exercício 3 %%%%%%%%%%%%%%%%%%%
% função SlidingWindow
function [DataWin, inipoint, endpoint, center] = SlidingWindow(Data, NwindowPoints, shift)
% Argumentos da função SlidingWindow
% Entradas:
%   Data          - vetor com o sinal no tempo
%   NwindowPoints - número de amostras de cada janela
%   shift         - deslocamento entre o início de uma janela e o da próxima
%
% Saídas:
%   DataWin  - matriz com NwindowPoints linhas e uma coluna por janela
%   inipoint - índice de início de cada janela
%   endpoint - índice de fim de cada janela
%   center   - índice do centro de cada janela

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
    idx = inipoint(k):min(endpoint(k), L);              % índices do sinal
    DataWin(1:length(idx), k) = Data(idx);              % restante fica zeros
    center(k) = ceil(median(inipoint(k):endpoint(k)));
end
end

%%%%%%%%%%%%%%%%%%% Exercício 4 %%%%%%%%%%%%%%%%%%%
function [Pxx, FreqVect] = WelchPeriodogram(Sinal, Fs, NFFT, NwindowPoints, shift, tipoJanela)
% WelchPeriodogram
%   [Pxx, FreqVect] = WelchPeriodogram(Sinal, Fs, NFFT, NwindowPoints, shift, tipoJanela)
%   tipoJanela - 'hamming' ou 'retangular'

DataWin = SlidingWindow(Sinal, NwindowPoints, shift);

if strcmpi(tipoJanela, 'hamming')
    w = hamming(NwindowPoints);         % janela de Hamming
else
    w = ones(NwindowPoints, 1);
end

Pxx = 0;
for k = 1:size(DataWin, 2)
    [Ymag, fase, FreqVect] = EspectAnalysis(DataWin(:, k).*w, Fs, NFFT, 0);
    Pxx = Pxx + Ymag.^2;
end
Pxx = Pxx/(size(DataWin, 2));         % média dos periodogramas das janelas
end

% Periodograma da realização 1 (POz)
NwindowPoints = 512;    % tamanho da janela (enunciado do exercicio)
shift = 256;            % deslocamento - sobreposição de 50 % (enunciado)

% Janela de Hamming
[PxxHam, fHam] = WelchPeriodogram(sinal1, fs, NFFT, NwindowPoints, shift, 'hamming');

% Janela retangular
[PxxRet, fRet] = WelchPeriodogram(sinal1, fs, NFFT, NwindowPoints, shift, 'retangular');

figure;

subplot(3,1,1);
plot(f1, Ymag1, 'k');
xlim([5 30]);
xlabel('Frequência (Hz)');
ylabel('|Y(f)|');
title('Realização 1 - Espectro de magnitude');
grid on;

subplot(3,1,2);
plot(fHam, PxxHam, 'k');
xlim([5 30]);
xlabel('Frequência (Hz)');
ylabel('Potência');
title('Realização 1 - Periodograma de Welch - Janela de Hamming');
grid on;

subplot(3,1,3);
plot(fRet, PxxRet, 'k');
xlim([5 30]);
xlabel('Frequência (Hz)');
ylabel('Potência');
title('Realização 1 - Periodograma de Welch - Janela retangular');
grid on;

% Pico entre 9 e 14 Hz em cada caso
faixaW = find(fHam >= 9 & fHam <= 14);
[valorHam, posHam] = max(PxxHam(faixaW));
[valorRet, posRet] = max(PxxRet(faixaW));
disp(['Realização 1: estímulo focado = ' num2str(frequencies(label(1))) ' Hz']);
disp(['Welch Hamming: pico em ' num2str(fHam(faixaW(posHam))) ' Hz']);
disp(['Welch retangular: pico em ' num2str(fRet(faixaW(posRet))) ' Hz']);
disp(['-----------------------------------------------------']);

% Resposta Final da 4
% Em relação ao espectro de magnitude, o periodograma de Welch apresenta um sinal bem mais suave,
% permitindo uma melhor observação do sinal.
% A Janela de Hamming produz um espectro mais suave, com picos mais largos e arredondados e menos oscilações entre eles, 
% já a janela retangular gera picos mais estreitos e definidos, mas o espectro fica mais irregular, com mais ondulações. 

% Na minha interpretação, a janela de retangular parece apresentar uma melhor "resolução", pois fica mais fácil identificar 
% os picos quando comparado com a janela de Hamming. 

%%%%%%%%%%%%%%%%%%% Exercício 5 %%%%%%%%%%%%%%%%%%%
% Comparação com a função pwelch do MATLAB (Signal Processing Toolbox)
[pxx, fw] = pwelch(sinal1, NwindowPoints, shift, NFFT, fs);

figure;

subplot(2,1,1);
plot(fHam, PxxHam, 'k');
xlim([5 30]);
xlabel('Frequência (Hz)');
ylabel('Potência');
title('Realização 1 - Welch implementado');
grid on;

subplot(2,1,2);
plot(fw, pxx, 'k');
xlim([5 30]);
xlabel('Frequência (Hz)');
ylabel('PSD');
title('Realização 1 - Função pwelch do MATLAB');
grid on;

% Resposta Final da 5
% As duas curvas apresentam o mesmo formato, independentemente da função utilizada.

%%%%%%%%%%%%%%%%%%% Exercício 6 %%%%%%%%%%%%%%%%%%%
% Classificação de todas as realizações pelo canal POz
freqEstimulo = [10 11 12 13];       % frequências dos estímulos
Ntrials = size(stimuliData, 3);     % número de realizações (48)
labelhat = zeros(Ntrials, 1);       % classes estimadas

for trial = 1:Ntrials
    sinal = stimuliData(:, canal, trial);
    [Pxx, f] = WelchPeriodogram(sinal, fs, NFFT, NwindowPoints, shift, 'hamming');

    % Potência no ponto mais próximo de cada frequência de estímulo
    potencia = zeros(1, 4);
    for j = 1:4
        [diferenca, posicao] = min(abs(f - freqEstimulo(j)));
        potencia(j) = Pxx(posicao);
    end

    % A classe é a frequência com maior potência
    [maiorPotencia, classe] = max(potencia);
    labelhat(trial) = classe;
end

taxaacerto = 1 - (sum(labelhat ~= label(:))/length(label));
disp(['Taxa de acerto (POz, Hamming) = ' num2str(taxaacerto*100) ' %']);
disp(['-----------------------------------------------------']);

% Resposta Final da 6
% Com o canal POz, janela de Hamming, NwindowPoints = 512, shift = 256 e NFFT = 2048,
% a taxa de acerto foi de 47,9 %.

%%%%%%%%%%%%%%%%%%% Exercício 7 %%%%%%%%%%%%%%%%%%%
% Exercício 7 - taxa de acerto de um canal, usando todas as realizações
function taxa = TaxaAcerto(stimuliData, label, canal, fs, NFFT, NwindowPoints, shift, tipoJanela)

freqEstimulo = [10 11 12 13];
Ntrials = size(stimuliData, 3);
labelhat = zeros(Ntrials, 1);

for trial = 1:Ntrials
    sinal = stimuliData(:, canal, trial);
    [Pxx, f] = WelchPeriodogram(sinal, fs, NFFT, NwindowPoints, shift, tipoJanela);

    potencia = zeros(1, 4);
    for j = 1:4
        [diferenca, posicao] = min(abs(f - freqEstimulo(j)));
        potencia(j) = Pxx(posicao);
    end

    [maiorPotencia, classe] = max(potencia);
    labelhat(trial) = classe;
end

taxa = 1 - (sum(labelhat ~= label(:))/length(label));
end

% parametros 
Ncanais = size(stimuliData, 2);     % número de canais (16)

% Parte a) e b): taxa de acerto de cada canal com janela de Hamming e retangular
taxaHam = zeros(Ncanais, 1);
taxaRet = zeros(Ncanais, 1);

for c = 1:Ncanais
    taxaHam(c) = TaxaAcerto(stimuliData, label, c, fs, NFFT, NwindowPoints, shift, 'hamming');
    taxaRet(c) = TaxaAcerto(stimuliData, label, c, fs, NFFT, NwindowPoints, shift, 'retangular');
end

disp('Taxa de acerto por canal (NFFT = 2048, janela = 512, shift = 256):');
disp('Canal   Hamming(%)   Retangular(%)');
for c = 1:Ncanais
    disp([num2str(c) '       ' num2str(taxaHam(c)*100) '        ' num2str(taxaRet(c)*100)]);
end

[picoHam, canalHam] = max(taxaHam);
[picoRet, canalRet] = max(taxaRet);
disp(['Melhor canal com Hamming: ' num2str(canalHam) ' (' num2str(picoHam*100) ' %)']);
disp(['Melhor canal com retangular: ' num2str(canalRet) ' (' num2str(picoRet*100) ' %)']);
disp(['-----------------------------------------------------']);

figure;
bar([taxaHam taxaRet]*100);
xlabel('Canal');
ylabel('Taxa de acerto (%)');
legend('Hamming', 'Retangular');
title('Taxa de acerto por canal');
grid on;

% Variação do NFFT, com janela = NFFT e sobreposição de 50 % (shift = NFFT/2)
listaNFFT = [256 512 1024 2048];
picoHamNFFT = zeros(1, 4);
picoRetNFFT = zeros(1, 4);

for b = 1:4
    NFFTb = listaNFFT(b);
    taxasHam = zeros(Ncanais, 1);
    taxasRet = zeros(Ncanais, 1);
    for c = 1:Ncanais
        taxasHam(c) = TaxaAcerto(stimuliData, label, c, fs, NFFTb, NFFTb, NFFTb/2, 'hamming');
        taxasRet(c) = TaxaAcerto(stimuliData, label, c, fs, NFFTb, NFFTb, NFFTb/2, 'retangular');
    end
    picoHamNFFT(b) = max(taxasHam);
    picoRetNFFT(b) = max(taxasRet);

    [janelas, ini, fim, centro] = SlidingWindow(stimuliData(:, 1, 1), NFFTb, NFFTb/2);
    disp(['NFFT = ' num2str(NFFTb) ' | janelas = ' num2str(size(janelas, 2)) ...
          ' | resolução = ' num2str(fs/NFFTb) ' Hz' ...
          ' | pico Hamming = ' num2str(picoHamNFFT(b)*100) ' %' ...
          ' | pico retangular = ' num2str(picoRetNFFT(b)*100) ' %']);
end

figure;
plot(listaNFFT, picoHamNFFT*100, '-ok'); hold on;
plot(listaNFFT, picoRetNFFT*100, '--sk'); hold off;
xticks(listaNFFT);
xlabel('NFFT (= tamanho da janela)');
ylabel('Taxa de acerto de pico (%)');
legend('Hamming', 'Retangular');
title('Desempenho de pico x NFFT (sobreposição de 50 %)');
grid on;

% Resposta Final da 7
% a) O canal 1 leva ao melhor desempenho (72,9 % com janela retangular e 66,6 % com
%    Hamming).
%
% b) A janela retangular leva ao melhor desempenho de pico (72,9 % contra 66,6 % da
%    Hamming).
%
% Desempenho de pico (melhor canal, sempre o canal 1) com sobreposição de 50 %:
%    NFFT = 256  -> 11 janelas, resolução de 1 Hz     -> Hamming 66,67 % | Retangular 70,83 %
%    NFFT = 512  -> 5 janelas,  resolução de 0,5 Hz   -> Hamming 66,67 % | Retangular 72,92 %
%    NFFT = 1024 -> 2 janelas,  resolução de 0,25 Hz  -> Hamming 75,00 % | Retangular 77,08 %
%    NFFT = 2048 -> 1 janela (com zeros), 0,125 Hz    -> Hamming 79,17 % | Retangular 75,00 %

% Existe um compromisso: janelas curtas dão mais janelas para fazer a média (menor variância,
% espectro mais estável), mas pior resolução espectral. Janelas longas dão melhor resolução,
% mas poucas janelas para a média (espectro mais ruidoso). Neste conjunto de dados, o
% desempenho melhorou com janelas maiores: a resolução espectral foi mais importante que o
% número de janelas, pois os estímulos estão separados por só 1 Hz. O melhor resultado foi
% 79,17 % (NFFT = 2048, Hamming), quando o sinal inteiro vira uma única janela. Nesse caso
% extremo não há média entre janelas, e a janela de Hamming passa a ser vantajosa porque
% reduz o vazamento espectral do sinal inteiro.


