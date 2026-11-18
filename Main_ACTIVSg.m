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
% Estudo de contingencias no sistema sintetico ACTIVSg2000 (Texas, 2000
% barras) com os comutadores de tap sob carga (LTCs) dos transformadores
% de rede, montados por Caso_ACTIVSg:
%   - horas 5368 (carga pesada) e 8504 (carga leve) do caso horario;
%   - perda de cada uma das 10 linhas de 500 kV mais carregadas da hora
%     (ramos com as duas barras em 500 kV; nao ha transformadores 500/500);
%   - os LTCs restabelecem a tensao das barras de baixa tensao com taps
%     discretos (32 posicoes em +-10%), banda morta e limites de reativo.
% HELM e NR usam o mesmo laco de eventos (HELMT_Resolver) e devem chegar
% aos mesmos eventos. Para medir o efeito dos LTCs, cada contingencia e
% resolvida tambem pelo HELM com os taps fixos no valor de base; anota-se
% o desvio |V| - Vesp nas barras controladas, sem e com os LTCs.
% Grava resultados/activsg.dat e, para a contingencia de maior queda de
% tensao sem os LTCs, os desvios de todas as barras controladas em
% resultados/activsg_perfil.dat.


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
disp('.    TRANSFORMADORES REGULADORES NO HELM: ACTIVSg2000 COM LTCs     .')
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


%% b) Limites e Degraus dos LTCs
lim.amin  = 0.90;
lim.amax  = 1.10;
lim.fmin  = 0;
lim.fmax  = 0;
lim.astep = 0.00625;
lim.fstep = 0;
lim.banda = 0.01;
lim.nmov  = 3;


%% c) Horas e Numero de Contingencias
horas  = [5368 8504];
nomesH = {'pesada', 'leve'};
ncont  = 10;
res    = zeros(numel(horas)*ncont, 16);
e0     = struct('ramo', {}, 'modo', {}, 'barra', {}, 'Vesp', {}, 'Psp', {});
pior   = struct('dV', 0, 'hora', 0, 'de', 0, 'para', 0, 'sem', [], 'com', []);
lin    = 0;



%% ________________________________________________________________________
%%
%%                         3. CONTINGENCIAS
%% ________________________________________________________________________
for ih = 1:numel(horas)

    %% a) Caso da Hora e LTCs
    [D, c, inf0] = Caso_ACTIVSg(horas(ih));
    V0 = inf0.V0;
    fprintf('\nHora %d (carga %s): %d LTCs em %d transformadores, %d barras PV passam a PQ\n', ...
        horas(ih), nomesH{ih}, numel(c), inf0.nRamos, inf0.nPVfixo);


    %% b) Linhas de 500 kV Mais Carregadas no Caso Base
    br   = D.branch;
    lin5 = find(D.kV(br(:, 1)) >= 500 & D.kV(br(:, 2)) >= 500);
    P5   = abs(FluxoT(D, lin5, V0, ones(numel(lin5), 1), zeros(numel(lin5), 1)))*D.baseMVA;
    [~, ord] = sort(P5, 'descend');

    % um circuito por corredor: circuitos paralelos iguais dao o mesmo caso
    par = sort(br(lin5(ord), 1:2), 2);
    [~, iu] = unique(par, 'rows', 'first');
    ord = ord(sort(iu));
    lst = lin5(ord(1:ncont));

    % tap base de cada grupo de LTC (primeiro transformador do grupo)
    a0 = zeros(numel(c), 1);
    for i = 1:numel(c)
        a0(i) = D.branch(c(i).ramo(1), 9);
    end

    fprintf('%-5s %-11s %7s %7s %6s %6s %6s %6s %9s %7s %7s\n', 'Cont', 'Linha', 'P (MW)', 'iguais', 'Q', 'taps', 'mov', 'passos', 'dif V(%)', 'tHELM', 'tNR');


    %% c) Abertura de Cada Linha
    for j = 1:ncont
        Dc = D;
        l  = lst(j);

        % a linha sai de servico com impedancia muito alta
        Dc.branch(l, 3) = 0;
        Dc.branch(l, 4) = 1e6;
        Dc.branch(l, 5) = 0;


        %% c1) NR e HELM com o Mesmo Laco de Eventos
        tic
        [Vn, an, fn, iN] = HELMT_Resolver(Dc, c, lim, 'NR', xi, nit, true);
        tN = toc;

        tic
        [Vh, ah, fh, iH] = HELMT_Resolver(Dc, c, lim, 'HELM', xi, nmax, true);
        tH = toc;


        %% c2) Conferencia
        igual = iN.conv && iH.conv && isequal(iN.satQ, iH.satQ) ...
            && isequal(iN.satA, iH.satA) && iN.passos == iH.passos ...
            && max(abs(ah - an)) < 1e-6;
        difV = max(abs(abs(Vh) - abs(Vn))./abs(Vn))*100;
        nQ   = sum(iH.satQ ~= 0);
        amov = sum(abs(ah(:) - a0) > 1e-9);
        nmov = sum(iH.nmov);


        %% c3) Desvio das Tensoes Controladas, sem e com os LTCs
        % sem LTC: os taps ficam no valor de base (sem controle)
        [Vs, as, fs, iS] = HELMT_Resolver(Dc, e0, lim, 'HELM', xi, nmax, true);
        bc   = [c.barra];
        ve   = [c.Vesp];
        dSem = abs(Vs(bc)) - ve(:);
        dCom = abs(Vh(bc)) - ve(:);

        % contingencia de maior queda sem os LTCs
        if min(dSem) < pior.dV
            pior.dV   = min(dSem);
            pior.hora = horas(ih);
            pior.de   = D.num(br(l, 1));
            pior.para = D.num(br(l, 2));
            pior.sem  = dSem;
            pior.com  = dCom;
        end


        %% c4) Resultados
        lin = lin + 1;
        res(lin, :) = [horas(ih), j, D.num(br(l, 1)), D.num(br(l, 2)), P5(ord(j)), igual, nQ, amov, nmov, ...
                       iH.passos, difV, min(abs(Vh)), min(dSem), max(dSem), min(dCom), max(dCom)];
        fprintf('%-5d %5d-%-5d %7.0f %7d %6d %6d %6d %6d %9.1e %7.2f %7.2f\n', ...
            j, D.num(br(l, 1)), D.num(br(l, 2)), P5(ord(j)), igual, nQ, amov, nmov, iH.passos, difV, tH, tN);
    end
end



%% ________________________________________________________________________
%%
%%                            4. GRAVACAO
%% ________________________________________________________________________

%% a) Resumo das Contingencias
fid = fopen(fullfile('resultados', 'activsg.dat'), 'w');
% os tempos ficam so na tela (dependem do computador)
fprintf(fid, 'hora cont de para P_MW iguais nQ taps_mov mov_disc passos difV_pct Vmin dVmin_sem dVmax_sem dVmin_com dVmax_com\n');
fprintf(fid, '%d %d %d %d %.1f %d %d %d %d %d %.3e %.5f %.5f %.5f %.5f %.5f\n', res.');
fclose(fid);


%% b) Perfil da Contingencia de Maior Queda (desvios em ordem crescente)
fprintf('\nMaior queda sem os LTCs: hora %d, corredor %d-%d (%.4f p.u.)\n', pior.hora, pior.de, pior.para, pior.dV);
fid = fopen(fullfile('resultados', 'activsg_perfil.dat'), 'w');
fprintf(fid, 'i dV_sem dV_com\n');
ds = sort(pior.sem);
dc = sort(pior.com);

for i = 1:numel(ds)
    fprintf(fid, '%d %.6f %.6f\n', i, ds(i), dc(i));
end

fclose(fid);
