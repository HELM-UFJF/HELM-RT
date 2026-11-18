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
% Este programa repete os estudos de 24 h com os limites operacionais:
%   - reativo dos geradores (barras PV);
%   - modulo do tap (lim.amin, lim.amax) e fase (lim.fmin, lim.fmax);
%   - taps e fases continuos e, em seguida, discretos (degrau lim.astep e
%     lim.fstep, banda morta de tensao lim.banda).
% O mesmo laco de eventos (HELMT_Resolver) e usado pelo HELM e pelo NR, e
% cada ponto e conferido: os dois metodos devem chegar ao mesmo conjunto
% de barras PV no limite, de taps e fases saturados e de passagens.
% A curva de carga e reduzida para pico de 1,0 vez a carga base.
% Para cada sistema grava resultados/limites_<sistema>.dat.


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
disp('.    TRANSFORMADORES REGULADORES NO HELM: LIMITES OPERACIONAIS     .')
disp('.                                                                  .')
disp('....................................................................')



%% ________________________________________________________________________
%%
%%                       2. DADOS E PARAMETROS
%% ________________________________________________________________________

%% a) Parametros dos Metodos
xi   = 1e-10;
nmax = 300;
nit  = 30;


%% b) Limites e Degraus
% amin, amax : faixa do tap; fmin, fmax : faixa da fase (rad)
% astep      : degrau do tap (0,625%, 32 posicoes em +-10%)
% fstep      : degrau da fase (1 grau); banda: banda morta de tensao (p.u.)
% nmov       : movimentos maximos de um degrau por transformador
lim.amin  = 0.95;
lim.amax  = 1.05;
lim.fmin  = -3*pi/180;
lim.fmax  = 3*pi/180;
lim.astep = 0.00625;
lim.fstep = pi/180;
lim.banda = 0.01;
lim.nmov  = 3;


%% c) Sistemas, Modos e Escala da Curva de Carga
sistemas = [30 57 118];
modos    = {'V', 'P', 'VP'};
nomes    = {'V', 'F', 'VF'};
casos    = {'I', 'II'};
tipos    = {'continuos', 'discretos'};
escala   = 1/1.4;



%% ________________________________________________________________________
%%
%%                        3. SIMULACOES DE 24 H
%% ________________________________________________________________________
for sis = sistemas

    %% a) Configuracao do Sistema
    cfg = Config_Sistema(sis);
    nh  = numel(cfg.carga);
    res = zeros(2*2*3*nh, 19);
    dV  = zeros(2*2*3*nh, 1);
    lin = 0;

    fprintf('\nSistema %d barras\n', sis);
    fprintf('%-9s %-5s %-4s %7s %7s %7s %7s %8s %9s %8s\n', 'Taps', 'Caso', 'Modo', 'iguais', 'h Q', 'h tap', 'h fase', 'passos', 'dif V(%)', 'HELM/NR');


    %% b) Continuo e Discreto, Casos de Vento e Modos
    for disc = [false true]
        for cs = 1:2
            for q = 1:3
                tH  = 0;
                tN  = 0;
                ini = lin + 1;

                for h = 1:nh

                    %% b1) Caso e Controles da Hora
                    [c, D] = Monta_Controles(cfg, escala*cfg.carga(h), cs, modos{q});


                    %% b2) NR e HELM com o Mesmo Laco de Eventos
                    tic
                    [Vn, an, fn, iN] = HELMT_Resolver(D, c, lim, 'NR', xi, nit, disc);
                    tN = tN + toc;

                    tic
                    [Vh, ah, fh, iH] = HELMT_Resolver(D, c, lim, 'HELM', xi, nmax, disc);
                    tH = tH + toc;


                    %% b3) Conferencia e Resultados do Ponto
                    % iguais: mesmos eventos, mesmas passagens e mesmos taps
                    igual = iN.conv && iH.conv && isequal(iN.satQ, iH.satQ) ...
                        && isequal(iN.satA, iH.satA) && isequal(iN.satF, iH.satF) ...
                        && iN.passos == iH.passos && max(abs(ah - an)) < 1e-6;
                    P1   = FluxoT(D, c(1).ramo, Vh, ah(1), fh(1))*D.baseMVA;
                    P2   = FluxoT(D, c(2).ramo, Vh, ah(2), fh(2))*D.baseMVA;
                    difV = max(abs(abs(Vh) - abs(Vn))./abs(Vn))*100;
                    lin  = lin + 1;

                    % [disc caso modo hora nQ satA1 satA2 satF1 satF2 a1 a2 f1 f2 V1 V2 P1 P2 passos iguais]
                    res(lin, :) = [disc, cs, q, h, sum(iH.satQ ~= 0), iH.satA(1), iH.satA(2), ...
                                   iH.satF(1), iH.satF(2), ah(1), ah(2), fh(1)*180/pi, fh(2)*180/pi, ...
                                   abs(Vh(c(1).barra)), abs(Vh(c(2).barra)), P1, P2, iH.passos, igual];
                    dV(lin) = difV;
                end


                %% b4) Resumo
                r = res(ini:lin, :);
                fprintf('%-9s %-5s %-4s %4d/%d %7d %7d %7d %8d %9.1e %8.2f\n', ...
                    tipos{disc + 1}, casos{cs}, nomes{q}, sum(r(:, 19)), nh, ...
                    sum(r(:, 5) > 0), sum(any(r(:, 6:7) ~= 0, 2)), sum(any(r(:, 8:9) ~= 0, 2)), ...
                    max(r(:, 18)), max(dV(ini:lin)), tH/tN);
            end
        end
    end


    %% c) Gravacao
    fid = fopen(fullfile('resultados', sprintf('limites_%d.dat', sis)), 'w');
    fprintf(fid, 'disc caso modo hora nQ satA1 satA2 satF1 satF2 a1 a2 f1 f2 V1 V2 P1 P2 passos iguais\n');
    fprintf(fid, '%d %d %d %d %d %d %d %d %d %.6f %.6f %.4f %.4f %.6f %.6f %.4f %.4f %d %d\n', res.');
    fclose(fid);
end
