function hpf = lerHPF(arquivo)
% LERHPF  Lê arquivos .hpf (formato HPF da Data Translation, usado pelo
% Delsys EMGworks Acquisition e pelo QuickDAQ) diretamente no MATLAB.
%
% Uso:
%   hpf = lerHPF('Sub-C005_tempo-0_tarefa-cvim_braço-nt_medida-torque_Plot_and_Store_Rep_1.0.hpf');
%   hpf.nomes          % nomes dos canais (quando presentes no XML)
%   hpf.fs             % frequências de amostragem encontradas no XML
%   hpf.dados{1}       % sinal do canal 1 (vetor coluna)
%   disp(hpf.xmlCanais{1})   % metadados completos dos canais
%
% Estrutura do arquivo: uma sequência de "chunks", cada um começando com
% int64 chunkID + int64 chunkSize (little-endian):
%   0x1000 = cabeçalho (XML)       0x2000 = informações dos canais (XML)
%   0x3000 = dados                 0x4000-0x6000 = eventos/índice (ignorados)
% As amostras são lidas como float32 (padrão do EMGworks). Se o sinal sair
% sem sentido, veja o campo de tipo de dado em hpf.xmlCanais e troque
% PRECISAO abaixo (ex.: 'int32' ou 'double').

PRECISAO = 'float32';
bytesAmostra = 4;   % 4 para float32/int32, 8 para double

fid = fopen(arquivo, 'r', 'ieee-le');
if fid < 0
    error('Não foi possível abrir o arquivo: %s', arquivo);
end
limpa = onCleanup(@() fclose(fid));

fseek(fid, 0, 'eof'); tamanho = ftell(fid); fseek(fid, 0, 'bof');

hpf = struct('arquivo', arquivo, 'xmlCabecalho', '', 'xmlCanais', {{}}, ...
             'nCanais', 0, 'nomes', {{}}, 'fs', [], 'dados', {{}});
pedacos = {};   % pedacos{canal} = lista de blocos de amostras
pos = 0;

while pos + 16 <= tamanho
    fseek(fid, pos, 'bof');
    id = fread(fid, 1, 'int64');
    sz = fread(fid, 1, 'int64');
    if isempty(sz) || sz <= 0 || pos + sz > tamanho
        break
    end

    switch id
        case 4096   % 0x1000 - cabeçalho
            fread(fid, 1, 'int32');   % creatorID
            fread(fid, 1, 'int64');   % versão do arquivo
            fread(fid, 1, 'int64');   % offset do chunk de índice
            hpf.xmlCabecalho = lerTexto(fid, pos + sz);

        case 8192   % 0x2000 - informações dos canais
            fread(fid, 1, 'int32');   % groupID
            hpf.nCanais = hpf.nCanais + fread(fid, 1, 'int32');
            hpf.xmlCanais{end+1} = lerTexto(fid, pos + sz);

        case 12288  % 0x3000 - dados
            fread(fid, 1, 'int32');   % groupID
            fread(fid, 1, 'int64');   % índice da primeira amostra
            n = fread(fid, 1, 'int32');           % nº de canais neste chunk
            info = fread(fid, [2 n], 'int32');    % [offset; nº de bytes]
            for k = 1:n
                fseek(fid, pos + info(1, k), 'bof');
                v = fread(fid, info(2, k) / bytesAmostra, PRECISAO);
                if numel(pedacos) < k, pedacos{k} = {}; end
                pedacos{k}{end+1} = v;
            end
    end
    pos = pos + sz;
end

hpf.dados = cellfun(@(c) vertcat(c{:}), pedacos, 'UniformOutput', false);

% Metadados úteis extraídos do XML (nomes de tags variam entre versões)
xml = strjoin(hpf.xmlCanais, newline);
t = regexp(xml, '<Name>([^<]*)</Name>', 'tokens');
hpf.nomes = cellfun(@(x) x{1}, t, 'UniformOutput', false);
t = regexp(xml, '<\w*(?:Frequency|Rate)\w*>([\d.,Ee+-]+)<', 'tokens');
hpf.fs = cellfun(@(x) str2double(strrep(x{1}, ',', '.')), t);

fprintf('%d canal(is) de dados lido(s):\n', numel(hpf.dados));
for k = 1:numel(hpf.dados)
    nome = '';
    if k <= numel(hpf.nomes), nome = hpf.nomes{k}; end
    fprintf('  canal %d %-25s %d amostras\n', k, nome, numel(hpf.dados{k}));
end
end

function s = lerTexto(fid, fim)
% Lê o restante do chunk como texto (ASCII ou UTF-16) e isola o XML.
b = fread(fid, fim - ftell(fid), '*uint8')';
if numel(b) > 1 && mean(b(2:2:end) == 0) > 0.4      % parece UTF-16LE
    b = b(1:2*floor(numel(b)/2));
    s = char(typecast(b, 'uint16'));
else
    s = char(b);
end
s(s == char(0)) = [];
i1 = find(s == '<', 1, 'first'); i2 = find(s == '>', 1, 'last');
if ~isempty(i1) && ~isempty(i2), s = s(i1:i2); end
end
