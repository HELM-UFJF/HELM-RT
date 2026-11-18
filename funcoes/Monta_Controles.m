function [c, D] = Monta_Controles(cfg, FC, cs, modo)
%MONTA_CONTROLES  Caso e controles de uma hora do estudo de 24 h.
%
%   [c, D] = Monta_Controles(cfg, FC, cs, modo)
%
%   Entradas
%     cfg  : configuracao do sistema (Config_Sistema)
%     FC   : multiplicador de carga da hora
%     cs   : caso de vento (1 = minimo, 2 = maximo)
%     modo : 'V', 'P' ou 'VP'
%
%   O fluxo especificado e cfg.pctP % do fluxo do caso base da mesma hora
%   (transformadores no tap do caso e fase nula), resolvido pelo NR.
%
%   Saidas
%     c : controles dos Transformadores 1 e 2
%     D : estrutura de dados da hora (Prepara_dados)

    %% a) Caso da Hora com os Parques Eolicos
    eol = struct('barra', {}, 'pct', {});

    for f = 1:2
        eol(f).barra = cfg.fazenda(f);
        eol(f).pct   = cfg.vento(cs, f);
    end

    D = Prepara_dados(cfg.SEP, FC, eol);


    %% b) Fluxo Base dos Transformadores (sem controle)
    e0 = struct('ramo', {}, 'modo', {}, 'barra', {}, 'Vesp', {}, 'Psp', {});
    V0 = NRT_Nucleo(D, e0, 1e-10, 30);


    %% c) Controles
    c = e0;

    for j = 1:2
        l          = cfg.ramo(j);
        c(j).ramo  = l;
        c(j).modo  = modo;
        c(j).barra = cfg.barra(j);
        c(j).Vesp  = cfg.Vesp(j);
        c(j).Psp   = cfg.pctP(j)/100*FluxoT(D, l, V0, D.branch(l, 9), D.fase0(l));
    end

end
