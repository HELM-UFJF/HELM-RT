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
% Tempos de processamento do HELM (HELMT_Nucleo) e do NR (NRT_Nucleo) em
% uma solucao, sem o laco de eventos:
%   - IEEE 30, 57 e 118 barras na hora de pico (fator de carga 1,4), Caso I
%     de vento, modos V, F e VF;
%   - ACTIVSg2000 nas horas 5368 e 8504 com os LTCs (modo V).
% Cada caso e resolvido nrep vezes e e anotada a mediana. Antes das
% medicoes, cada caso e resolvido uma vez para aquecer o MATLAB.
% Grava resultados/benchmark_<versao>.dat.


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
disp('.      TRANSFORMADORES REGULADORES NO HELM: TEMPOS DE SOLUCAO      .')
disp('.                                                                  .')
disp('....................................................................')



%% ________________________________________________________________________
%%
%%                       2. DADOS E PARAMETROS
%% ________________________________________________________________________

%% a) Parametros dos Metodos
% mesmos criterios dos estudos: 1e-10 nos IEEE e no ACTIVSg2000
xi   = 1e-10;
nmax = 300;
nit  = 30;


%% b) Repeticoes
nrep  = 50;
nrepG = 10;


%% c) Casos
sistemas = [30 57 118];
modos    = {'V', 'P', 'VP'};
nomes    = {'V', 'F', 'VF'};
rot      = {};
tab      = [];



%% ________________________________________________________________________
%%
%%                          3. SISTEMAS IEEE
%% ________________________________________________________________________
for sis = sistemas
    cfg = Config_Sistema(sis);

    for q = 1:3

        %% a) Caso na Hora de Pico
        [c, D] = Monta_Controles(cfg, max(cfg.carga), 1, modos{q});


        %% b) Aquecimento
        [Vh, ah, fh, iH] = HELMT_Nucleo(D, c, xi, nmax);
        [Vn, an, fn, iN] = NRT_Nucleo(D, c, xi, nit);


        %% c) Medicoes
        tH = zeros(nrep, 1);
        tN = zeros(nrep, 1);

        for r = 1:nrep
            tic
            [Vh, ah, fh, iH] = HELMT_Nucleo(D, c, xi, nmax);
            tH(r) = toc;

            tic
            [Vn, an, fn, iN] = NRT_Nucleo(D, c, xi, nit);
            tN(r) = toc;
        end

        rot{end+1} = sprintf('IEEE %d %s', sis, nomes{q});
        tab(end+1, :) = [median(tH), median(tN), iH.n, iN.it];
    end
end



%% ________________________________________________________________________
%%
%%                           4. ACTIVSg2000
%% ________________________________________________________________________
for hora = [5368 8504]

    %% a) Caso com os LTCs
    [D, c] = Caso_ACTIVSg(hora);


    %% b) Aquecimento
    [Vh, ah, fh, iH] = HELMT_Nucleo(D, c, xi, nmax);
    [Vn, an, fn, iN] = NRT_Nucleo(D, c, xi, nit);


    %% c) Medicoes
    tH = zeros(nrepG, 1);
    tN = zeros(nrepG, 1);

    for r = 1:nrepG
        tic
        [Vh, ah, fh, iH] = HELMT_Nucleo(D, c, xi, nmax);
        tH(r) = toc;

        tic
        [Vn, an, fn, iN] = NRT_Nucleo(D, c, xi, nit);
        tN(r) = toc;
    end

    rot{end+1} = sprintf('ACTIVSg h%d', hora);
    tab(end+1, :) = [median(tH), median(tN), iH.n, iN.it];
end



%% ________________________________________________________________________
%%
%%                     5. RESULTADOS E GRAVACAO
%% ________________________________________________________________________

%% a) Tabela na Tela
fprintf('\n%s\n', version);
fprintf('%-16s %10s %10s %8s %7s %7s\n', 'Caso', 'HELM (s)', 'NR (s)', 'HELM/NR', 'termos', 'it NR');

for i = 1:numel(rot)
    fprintf('%-16s %10.4f %10.4f %8.2f %7d %7d\n', rot{i}, tab(i, 1), tab(i, 2), tab(i, 1)/tab(i, 2), tab(i, 3), tab(i, 4));
end


%% b) Arquivo
ver = regexprep(version, '[^0-9A-Za-z\.]', '_');
fid = fopen(fullfile('resultados', sprintf('benchmark_%s.dat', ver)), 'w');
fprintf(fid, '%s\n', version);
fprintf(fid, 'caso;HELM_s;NR_s;termos;itNR\n');

for i = 1:numel(rot)
    fprintf(fid, '%s;%.6f;%.6f;%d;%d\n', rot{i}, tab(i, 1), tab(i, 2), tab(i, 3), tab(i, 4));
end

fclose(fid);
