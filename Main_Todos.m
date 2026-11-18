%                                                                       JMJ
%
%       "Aquilo que Deus quer para nos e muito mais belo do que tudo
%              que poderiamos pedir com a nossa imaginacao."
%
%                                 Chiara Corbella Petrillo


%% ________________________________________________________________________
%%
%%                        1. INFORMACOES INICIAIS
%% ________________________________________________________________________

%% a) Comentarios Iniciais
% Este programa reproduz todos os resultados do artigo e confere cada
% arquivo gerado em resultados/ com o arquivo de mesmo nome em
% referencia/ (valores do artigo).
%   - Roda, em sequencia: Main_Modos, Main_Limites, Main_ACTIVSg e
%     Main_Convergencia;
%   - Main_Benchmark nao e rodado aqui, porque os tempos dependem do
%     computador; rode-o separadamente;
%   - Cada programa principal limpa as variaveis (clear all); por isso a
%     sequencia e escrita sem laco e a conferencia e feita no final, so a
%     partir dos arquivos gravados.
% Tempo total aproximado: alguns minutos.


%% b) Limpeza de Figuras e Variaveis
clc
clear all
close all
format short
warning off



%% ________________________________________________________________________
%%
%%                     2. EXECUCAO DOS PROGRAMAS PRINCIPAIS
%% ________________________________________________________________________

%% a) Modos de Controle sem Limites, IEEE 30, 57 e 118 Barras
Main_Modos


%% b) Limites Operacionais e Taps Discretos, IEEE 30, 57 e 118 Barras
Main_Limites


%% c) Contingencias no ACTIVSg2000 com LTCs
Main_ACTIVSg


%% d) Convergencia do Algoritmo Epsilon de Wynn
Main_Convergencia



%% ________________________________________________________________________
%%
%%                   3. CONFERENCIA COM OS VALORES DO ARTIGO
%% ________________________________________________________________________

%% a) Tolerancia e Lista de Arquivos de Referencia
% Um valor confere quando |gerado - referencia| <= tolAbs + tolRel*|referencia|.
clc
tolAbs = 1e-6;
tolRel = 1e-6;
lista  = dir(fullfile('referencia', '*.dat'));
nOk    = 0;
nFalha = 0;
nFalta = 0;

disp('....................................................................')
disp('.                                                                  .')
disp('.       TRANSFORMADORES REGULADORES: CONFERENCIA COM O ARTIGO      .')
disp('.                                                                  .')
disp('....................................................................')
fprintf('\n%-30s %12s   %s\n', 'Arquivo', 'max|dif|', 'Situacao');


%% b) Comparacao Arquivo a Arquivo
for i = 1:numel(lista)

    %% b1) Arquivo Gerado
    nome = lista(i).name;
    arqG = fullfile('resultados', nome);

    if ~exist(arqG, 'file')
        fprintf('%-30s %12s   %s\n', nome, '-', 'nao gerado');
        nFalta = nFalta + 1;
        continue
    end


    %% b2) Leitura dos Dois Arquivos (cabecalho na primeira linha)
    fid = fopen(fullfile('referencia', nome), 'r');
    cab = fgetl(fid);
    R   = textscan(fid, '%f');
    fclose(fid);

    fid = fopen(arqG, 'r');
    fgetl(fid);
    G   = textscan(fid, '%f');
    fclose(fid);

    R = R{1};
    G = G{1};


    %% b3) Diferenca (NaN nas mesmas posicoes conta como igual)
    if numel(R) ~= numel(G)
        fprintf('%-30s %12s   %s\n', nome, '-', 'tamanho diferente');
        nFalha = nFalha + 1;
        continue
    end

    ambosNaN = isnan(R) & isnan(G);
    dif      = abs(G - R);
    dif(ambosNaN) = 0;
    lim      = tolAbs + tolRel*abs(R);
    lim(ambosNaN) = 1;

    if all(dif <= lim)
        fprintf('%-30s %12.1e   %s\n', nome, max(dif), 'confere');
        nOk = nOk + 1;
    else
        fprintf('%-30s %12.1e   %s\n', nome, max(dif), 'DIVERGE');
        nFalha = nFalha + 1;
    end
end


%% c) Resumo
fprintf('\n%d arquivo(s) conferem, %d divergem, %d nao foram gerados.\n', nOk, nFalha, nFalta);

