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
% Convergencia do HELM com o algoritmo epsilon de Wynn: para cada caso,
% grava a variacao maxima da estimativa entre duas ordens consecutivas,
% epsilon_n = max |S_hat(n) - S_hat(n-1)|, sobre todas as series (tensoes,
% desvios de tap e fases).
%   - IEEE 30, 57 e 118 barras na hora de pico, Caso I, modo VF;
%   - ACTIVSg2000 nas horas 5368 e 8504 com os LTCs (modo V).
% Para comparacao, grava tambem o maior termo de cada ordem, que e a
% variacao da soma parcial simples (sem aceleracao).
% Grava resultados/convergencia.dat (duas colunas por caso, NaN apos a
% parada).


%% b) Limpeza de Figuras e Variaveis
clc
clear all
close all
format short
warning off


%% c) Pastas das Funcoes e dos Dados
addpath('funcoes', 'dados');

if ~exist('resultados', 'dir')
    mkdir('resultados');
end


%% d) Cabecalho
disp('....................................................................')
disp('.                                                                  .')
disp('.    TRANSFORMADORES REGULADORES NO HELM: CONVERGENCIA DO WYNN     .')
disp('.                                                                  .')
disp('....................................................................')



%% ________________________________________________________________________
%%
%%                       2. DADOS E PARAMETROS
%% ________________________________________________________________________

%% a) Parametros do HELM
xi   = 1e-10;
nmax = 300;


%% b) Casos
nomes = {'IEEE30_VF', 'IEEE57_VF', 'IEEE118_VF', 'ACTIVSg_h5368', 'ACTIVSg_h8504'};
H     = cell(1, 5);
G     = cell(1, 5);



%% ________________________________________________________________________
%%
%%                          3. SIMULACOES
%% ________________________________________________________________________

%% a) Sistemas IEEE na Hora de Pico
sistemas = [30 57 118];

for i = 1:3
    cfg    = Config_Sistema(sistemas(i));
    [c, D] = Monta_Controles(cfg, max(cfg.carga), 1, 'VP');
    [V, a, f, inf1] = HELMT_Nucleo(D, c, xi, nmax);
    H{i}   = inf1.hist;
    G{i}   = inf1.histS;
end


%% b) ACTIVSg2000 com os LTCs
horas = [5368 8504];

for i = 1:2
    [D, c] = Caso_ACTIVSg(horas(i));
    [V, a, f, inf1] = HELMT_Nucleo(D, c, xi, nmax);
    H{3+i} = inf1.hist;
    G{3+i} = inf1.histS;
end



%% ________________________________________________________________________
%%
%%                     4. RESULTADOS E GRAVACAO
%% ________________________________________________________________________

%% a) Matriz com uma Coluna por Caso
% colunas: [Wynn caso 1, soma parcial caso 1, Wynn caso 2, ...]
nl = max(cellfun(@numel, H));
M  = nan(nl, 10);

for i = 1:5
    M(1:numel(H{i}), 2*i-1) = H{i};
    M(1:numel(G{i}), 2*i)   = G{i};
end


%% b) Tela
fprintf('\n%-15s %8s %12s\n', 'Caso', 'ordens', 'epsilon_n');

for i = 1:5
    fprintf('%-15s %8d %12.2e\n', nomes{i}, numel(H{i}), H{i}(end));
end


%% c) Arquivo
fid = fopen(fullfile('resultados', 'convergencia.dat'), 'w');
fprintf(fid, 'n');

for i = 1:5
    fprintf(fid, ' %s_W %s_S', nomes{i}, nomes{i});
end

fprintf(fid, '\n');

for n = 1:nl
    fprintf(fid, '%d', n);
    fprintf(fid, ' %.6e', M(n, :));
    fprintf(fid, '\n');
end

fclose(fid);
