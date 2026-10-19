% Example 3.6 
% Determine and plot the power spectra of heart rate variability during both normal and meditative states.

% Solution 
% The PS can be obtained through the Fourier transform using the direct method given in Equation 3.34. 
% However, the heart rate data should first be converted into evenly sampled time data and this is a bit 
% tricky. The data set obtained by downloading from the PhysioNet database provides the heart rate at 
% unevenly spaced times, where the sample times are provided as a second vector. These interval data need 
% to be rearranged into evenly spaced time positions. This process, known as resampling, is done through 
% interpolation using MATLAB’s interp1 routine. This routine takes in the unevenly spaced x–y pairs as two 
% vectors along with a vector containing the desired evenly spaced x values. The routine then uses linear 
% interpolation (other options are possible) to approximate the y values that match the evenly spaced 
% x values. The details can be found in the MATLAB help file for interp1.


% In the program below, the uneven x–y pairs for the normal conditions are in vectors t_pre 
% and hr_pre, respectively, both in the MATLAB file Hr_pre. (For meditative conditions, the 
% vectors are named t_med and hr_med and are in file Hr_med). The evenly spaced time vector 
% is xi and the resampled heart rate values are in vector yi.
clear all;
close all;
clc;

addpath(genpath("funcoes_auxiliares_and_data_tests"));   % torna sig_noise visível

% Example 3.6
% Frequency analysis of heart rate data in the normal and meditative state
fs = 100;               % sample frequency (100 Hz)
ts = 1/fs;              % Sample interval
load Hr_pre;            % Load normal and meditative data

%
% Convert to evenly-spaced time data using interpolation; i.e., resampling
% First generate evenly space time vectors having one second
% intervals and extending over the time range of the data
%
xi = (ceil(t_pre(1)):ts:floor(t_pre(end)));     % Evenly-spaced time vector
yi = interp1(t_pre,hr_pre,xi');                 % Interpolate
yi = yi - mean(yi);                             % remove average
N2 = round(length(yi)/2);
f = (1:N2)*fs/(2*N2);                               % vector for plotting

% Now determine the Power spectrum
YI = abs((fft(yi)).^2); 

subplot(1,2,1);
plot(f,YI(2:N2 + 1),'k');
axis([0 .15 0 max(YI)*1.25]); 
title('HR - normal condition')
xlabel('Frequency (Hz)')
ylabel('Power Spectrum')
%
% Repeat for meditative data
%
fs = 100;               % sample frequency (100 Hz)
ts = 1/fs;              % Sample interval
load Hr_med;            % Load normal and meditative data

%
% Convert to evenly-spaced time data using interpolation; i.e., resampling
% First generate evenly space time vectors having one second
% intervals and extending over the time range of the data
%
xi = (ceil(t_med(1)):ts:floor(t_med(end)));     % Evenly-spaced time vector
yi = interp1(t_med,hr_med,xi');                 % Interpolate
yi = yi - mean(yi);                             % remove average
N2 = round(length(yi)/2);
f = (1:N2)*fs/(2*N2);                               % vector for plotting

% Now determine the Power spectrum
YI = abs((fft(yi)).^2); 

subplot(1,2,2);
plot(f,YI(2:N2 + 1),'k');
axis([0 .15 0 max(YI)*1.25]); 
title('HR - meditative')
xlabel('Frequency (Hz)')
ylabel('Power Spectrum')
%
% Análise
% Para converter os dados da frequência cardíaca em uma sequência de pontos uniformemente espaçados no tempo, 
% cria-se primeiramente um vetor temporal, x_i, que aumenta em incrementos de 0,01 s (1/f_s = 1/100) entre 
% os valores mínimo e máximo de tempo (devidamente arredondados) presentes nos dados originais. Utiliza-se 
% aqui uma frequência de reamostragem de 100 Hz, uma vez que esse valor é comum em estudos de variabilidade 
% da frequência cardíaca que empregam certos métodos não lineares; contudo, neste exemplo, uma ampla gama de 
% frequências de reamostragem produziria o mesmo resultado. 
% Os dados temporais uniformemente espaçados, y_i, foram gerados utilizando a rotina de interpolação 
% `interp1` do MATLAB. Este exemplo solicita o espectro de potência (PS, na sigla em inglês) da 
% variabilidade da frequência cardíaca — e não da frequência cardíaca em si —, ou seja, a variação da 
% frequência batimento a batimento. Para obter essa variação batimento a batimento, é necessário calcular a 
% diferença entre os valores das amostras, utilizando o operador ` yi = yi - mean(yi);` no MATLAB, antes de 
% determinar o espectro de potência.
% Após a interpolação e a remoção da frequência cardíaca média, o espectro de potência é determinado por meio 
% da transformada rápida de Fourier (`fft`), elevando-se, em seguida, o componente de magnitude ao quadrado. 
% Os gráficos de frequência estão limitados a um máximo de 0,15 Hz, visto que é nessa faixa que se concentra 
% a maior parte da energia espectral. Observe que o componente de corrente contínua (termo DC) não é 
% representado no gráfico.

% Resultados: 
% O espectro de potência da variabilidade da frequência cardíaca em condições normais apresenta valores 
% baixos, com um leve aumento em função da frequência. O espectro obtido durante a meditação exibe um pico 
% acentuado em torno de 0,12 Hz, o que indica a presença de um processo ressonante nessas frequências 
% — correspondendo a um intervalo de tempo de aproximadamente 8 segundos. 


%%%%%%%%%%%% EXPLICAÇÃO EXTRA %%%%%%%%%%%%
% As duas figuras mostram o espectro de potência da variabilidade da frequência cardíaca (VFC), ou seja, como as oscilações da 
% FC se distribuem entre as frequências, em repouso normal (esquerda) e durante meditação (direita). Como a FC varia lentamente, 
% as frequências envolvidas são muito baixas: 0,1 Hz corresponde a um ciclo a cada 10 segundos.
% Condição normal
% A potência está espalhada, sem um ritmo dominante bem definido. Ela se concentra nas frequências mais baixas, abaixo de ~0,05 Hz, 
% e vai diminuindo à medida que a frequência aumenta. O maior pico, em ~0,017 Hz, corresponde a uma oscilação com período de ~60 s, 
% mas ele é apenas um entre vários picos de altura parecida. Esse padrão é típico da VFC em repouso: a FC flutua de forma irregular, 
% influenciada por vários mecanismos ao mesmo tempo, como a termorregulação, o sistema renina-angiotensina e o barorreflexo, sem que 
% nenhum deles domine.
% Isso é coerente com a autocovariância da FC em repouso que você viu antes: correlação entre batimentos próximos que decai 
% suavemente, sem picos repetidos, ou seja, sem periodicidade forte.

% Meditação
% O padrão é completamente diferente. Quase toda a potência se concentra em um pico estreito e dominante em ~0,105 Hz, que corresponde
%  a um período de ~9,5 s. A FC passa a oscilar de forma regular, subindo e descendo a cada ~10 s. A potência nas frequências muito 
%  baixas praticamente desaparece.
% A explicação fisiológica mais aceita é a respiração lenta. Nas técnicas de meditação desses dados, a respiração desacelera para 
% cerca de 6 ciclos por minuto, ou seja, 0,1 Hz. A FC acompanha a respiração (aumenta na inspiração e diminui na expiração), um 
% fenômeno chamado arritmia sinusal respiratória. Por volta de 0,1 Hz, esse efeito entra em ressonância com o barorreflexo, que 
% tem frequência natural próxima disso. O resultado são oscilações da FC muito amplas e regulares, que aparecem como esse pico 
% concentrado.
% Cuidado ao comparar as escalas
% Os eixos verticais são diferentes: ×10⁹ na condição normal e ×10¹⁰ na meditação. O pico da meditação (~2,9 × 10¹⁰) é cerca de 5 
% vezes maior que o maior pico da condição normal (~6 × 10⁹). Ou seja, a meditação não só reorganiza a potência em uma única 
% frequência, como também aumenta muito a amplitude das oscilações.
% Essa comparação de valores absolutos tem uma ressalva. Como o espectro foi calculado sem normalização, a potência cresce com o 
% número de amostras. Se os dois registros tiverem durações diferentes, parte da diferença de escala vem só do comprimento dos dados. 
% Para uma comparação rigorosa, normalize dividindo por N, ou compare a distribuição relativa da potência entre as faixas.

% Relação com as bandas clássicas da VFC
% Na análise de VFC, o espectro costuma ser dividido em três bandas:
% •	VLF (< 0,04 Hz): onde está a maior parte da potência na condição normal;
% •	LF (0,04–0,15 Hz): onde está o pico dominante da meditação;
% •	HF (0,15–0,4 Hz): onde normalmente aparece a respiração espontânea, com frequências de ~0,2 a 0,3 Hz.
% Como o eixo foi limitado a 0,15 Hz, a banda HF ficou de fora. Em repouso normal, é provável que exista um pico respiratório nessa 
% faixa que não aparece no gráfico. Se ampliar o eixo para axis([0 0.4 ...]), você pode comparar também essa região. Na meditação, o 
% que acontece é que a respiração sai da HF e entra na LF, porque ficou mais lenta.
% Em resumo
% Em repouso, a FC tem flutuações lentas e irregulares, com potência distribuída entre várias frequências baixas. Na meditação, a 
% respiração lenta impõe uma oscilação regular de ~0,1 Hz, que concentra a potência em um único pico muito mais intenso. O exemplo 
% mostra como o espectro de potência revela mudanças na dinâmica de um sinal fisiológico que seriam difíceis de perceber olhando só 
% para a série temporal.
