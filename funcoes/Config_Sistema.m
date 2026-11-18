function cfg = Config_Sistema(sis)
%CONFIG_SISTEMA  Configuracao dos estudos de 24 h nos sistemas IEEE.
%
%   cfg = Config_Sistema(sis)      sis = 30, 57 ou 118
%
%   Campos de cfg
%     SEP     : caso no formato MATPOWER
%     ramo    : ramos dos Transformadores 1 e 2 (lado k = DE, com o tap)
%     barra   : barras controladas em tensao (PQ)
%     Vesp    : tensoes especificadas (p.u.)
%     pctP    : fluxo especificado em % do fluxo base (modos F e VF)
%     fazenda : barras dos Parques Eolicos A e B (potencia nominal igual a
%               carga ativa da barra)
%     carga   : multiplicador horario da carga (24 h)
%     vento   : % da potencia nominal dos parques nos Casos I e II
%               (cenarios representativos de minima e maxima geracao)

    %% a) Transformadores Controlados e Parques Eolicos
    switch sis
        case 30
            cfg.SEP     = case30();
            cfg.ramo    = [11 36];
            cfg.barra   = [9 27];
            cfg.Vesp    = [1.035 1.025];
            cfg.fazenda = [7 21];

        case 57
            cfg.SEP     = case57();
            cfg.ramo    = [71 46];
            cfg.barra   = [43 32];
            cfg.Vesp    = [1.030 0.980];
            cfg.fazenda = [15 29];

        case 118
            cfg.SEP     = case118();
            cfg.ramo    = [36 51];
            cfg.barra   = [17 37];
            cfg.Vesp    = [0.995 0.990];
            cfg.fazenda = [60 78];

        otherwise
            error('Config_Sistema: sistema %d nao configurado.', sis);
    end


    %% b) Fluxo Especificado (% do fluxo base), Curva de Carga e Vento
    cfg.pctP  = [120 80];
    cfg.carga = [0.70, 0.70, 0.70, 0.75, 0.80, 0.90, ...
                 1.00, 1.10, 1.20, 1.30, 1.35, 1.40, ...
                 1.35, 1.30, 1.25, 1.20, 1.15, 1.20, ...
                 1.25, 1.20, 1.10, 1.00, 0.90, 0.80];
    cfg.vento = [ 1.71,  3.36; ...
                 96.28, 98.47];

end
