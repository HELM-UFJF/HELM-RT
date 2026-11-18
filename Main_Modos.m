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
% Este programa resolve os tres modos de controle dos transformadores
% reguladores sem limites operacionais:
%   - V  : o modulo do tap regula a tensao de uma barra PQ;
%   - F  : a fase regula o fluxo ativo no lado k do transformador;
%   - VF : os dois controles ao mesmo tempo.
% Sistemas IEEE 30, 57 e 118, com dois transformadores controlados, dois
% parques eolicos e a curva de carga de 24 h, nos Casos I (vento minimo) e
% II (vento maximo). Cada ponto e resolvido pelo HELM (Wynn) e pelo NR
% completo, com a mesma tolerancia, e os dois sao comparados.
% Para cada sistema grava resultados/modos_<sistema>.dat.


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
disp('.      TRANSFORMADORES REGULADORES NO HELM: MODOS DE CONTROLE      .')
disp('.                                                                  .')
disp('....................................................................')



%% ________________________________________________________________________
%%
%%                       2. DADOS E PARAMETROS
%% ________________________________________________________________________

%% a) Parametros
% xi   : tolerancia do HELM (Wynn) e do NR
% nmax : ordem maxima do HELM
% nit  : iteracoes maximas do NR
xi   = 1e-10;
nmax = 300;
nit  = 30;


%% b) Sistemas e Modos
sistemas = [30 57 118];
modos    = {'V', 'P', 'VP'};
nomes    = {'V', 'F', 'VF'};
casos    = {'I', 'II'};



%% ________________________________________________________________________
%%
%%                        3. SIMULACOES DE 24 H
%% ________________________________________________________________________
for sis = sistemas

    %% a) Configuracao do Sistema
    cfg = Config_Sistema(sis);
    nh  = numel(cfg.carga);
    nb  = size(cfg.SEP.bus, 1);
    res = zeros(2*3*nh, 16);
    lin = 0;

    fprintf('\nSistema %d barras\n', sis);
    fprintf('%-6s %-4s %6s %10s %11s %11s %13s %13s %8s %8s\n', 'Caso', 'Modo', 'ok', 'termos', 'tap 1', 'tap 2', 'fase 1 (gr)', 'fase 2 (gr)', 'dif V(%)', 'HELM/NR');


    %% b) Casos de Vento, Modos e Horas
    for cs = 1:2
        for q = 1:3
            tH = 0;
            tN = 0;
            ini = lin + 1;

            for h = 1:nh

                %% b1) Caso e Controles da Hora
                [c, D] = Monta_Controles(cfg, cfg.carga(h), cs, modos{q});


                %% b2) NR e HELM
                tic
                [Vn, an, fn, iN] = NRT_Nucleo(D, c, xi, nit);
                tN = tN + toc;

                tic
                [Vh, ah, fh, iH] = HELMT_Nucleo(D, c, xi, nmax);
                tH = tH + toc;


                %% b3) Resultados do Ponto
                % [caso modo hora a1 a2 f1 f2 |Vc1| |Vc2| P1 P2 Vmin termos itNR difV ok]
                P1   = FluxoT(D, c(1).ramo, Vh, ah(1), fh(1))*D.baseMVA;
                P2   = FluxoT(D, c(2).ramo, Vh, ah(2), fh(2))*D.baseMVA;
                difV = max(abs(abs(Vh) - abs(Vn))./abs(Vn))*100;
                ok   = iN.conv && iH.conv;
                lin  = lin + 1;
                res(lin, :) = [cs, q, h, ah(1), ah(2), fh(1)*180/pi, fh(2)*180/pi, ...
                               abs(Vh(c(1).barra)), abs(Vh(c(2).barra)), P1, P2, ...
                               min(abs(Vh)), iH.n, iN.it, difV, ok];
            end


            %% b4) Resumo do Caso e do Modo
            r = res(ini:lin, :);
            fprintf('%-6s %-4s %3d/%d %4d a %3d %4.3f-%4.3f %4.3f-%4.3f %5.2f/%5.2f %5.2f/%5.2f %8.1e %8.2f\n', ...
                casos{cs}, nomes{q}, sum(r(:, 16)), nh, min(r(:, 13)), max(r(:, 13)), ...
                min(r(:, 4)), max(r(:, 4)), min(r(:, 5)), max(r(:, 5)), ...
                min(r(:, 6)), max(r(:, 6)), min(r(:, 7)), max(r(:, 7)), max(r(:, 15)), tH/tN);
        end
    end


    %% c) Gravacao
    fid = fopen(fullfile('resultados', sprintf('modos_%d.dat', sis)), 'w');
    fprintf(fid, 'caso modo hora a1 a2 f1 f2 V1 V2 P1 P2 Vmin termos itNR difV ok\n');
    fprintf(fid, '%d %d %d %.6f %.6f %.4f %.4f %.6f %.6f %.4f %.4f %.6f %d %d %.3e %d\n', res.');
    fclose(fid);
end
