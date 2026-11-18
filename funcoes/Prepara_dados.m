function D = Prepara_dados(SEP, FC, eol)
%PREPARA_DADOS  Organiza um caso no formato MATPOWER na estrutura D usada
%   pelo HELM e pelo NR, no fator de carregamento FC.
%
%   D = Prepara_dados(SEP, FC, eol)
%
%   Entradas
%     SEP : caso (SEP.baseMVA, SEP.bus, SEP.gen, SEP.branch); as barras podem
%           ter numeracao qualquer, que e convertida para 1..nb
%     FC  : fator de carregamento (multiplica as cargas ativa e reativa)
%     eol : (opcional) geracao eolica, vetor de estruturas com
%           .barra  numero da barra (deve ser PQ)
%           .pct    potencia gerada em % da demanda ativa nominal da barra
%           (injecao ativa com fator de potencia unitario)
%
%   Saida (campos de D)
%     baseMVA, nb, tipo (1 PQ, 2 PV, 3 V-theta), ref, pv, pq, Vesp, Sesp
%     (geracao - carga, p.u.), Qd, Qmin, Qmax (barras PV, p.u.), bus,
%     branch (barras renumeradas, ramos ligados, tap 0 -> 1), fase0 (fase
%     de cada ramo, rad), num (numeros originais), kV

    %% ____________________________________________________________________
    %%
    %%                     1. LEITURA E RENUMERACAO
    %% ____________________________________________________________________

    %% a) Dados de Entrada
    baseMVA = SEP.baseMVA;
    bus     = SEP.bus;
    gen     = SEP.gen;
    branch  = SEP.branch(:, 1:11);
    nb      = size(bus, 1);


    %% b) Ramos e Geradores Ligados
    branch = branch(branch(:, 11) > 0, :);
    gen    = gen(gen(:, 8) > 0, :);


    %% c) Numeracao Sequencial das Barras
    num = bus(:, 1);
    mp  = zeros(max(num), 1);
    mp(num) = 1:nb;
    bus(:, 1)    = (1:nb).';
    branch(:, 1) = mp(branch(:, 1));
    branch(:, 2) = mp(branch(:, 2));
    gen(:, 1)    = mp(gen(:, 1));


    %% d) Tipo de Barra
    tipo = bus(:, 2);
    ref  = find(tipo == 3);
    pv   = find(tipo == 2);
    pq   = find(tipo == 1);



    %% ____________________________________________________________________
    %%
    %%                   2. GERACAO, CARGA E LIMITES
    %% ____________________________________________________________________

    %% a) Geracao por Barra
    Pg   = accumarray(gen(:, 1), gen(:, 2), [nb, 1])/baseMVA;
    Qg   = accumarray(gen(:, 1), gen(:, 3), [nb, 1])/baseMVA;
    Qmx  = accumarray(gen(:, 1), gen(:, 4), [nb, 1])/baseMVA;
    Qmn  = accumarray(gen(:, 1), gen(:, 5), [nb, 1])/baseMVA;

    % Limites de reativo apenas nas barras PV
    Qmax = inf(nb, 1);
    Qmin = -inf(nb, 1);
    Qmax(pv) = Qmx(pv);
    Qmin(pv) = Qmn(pv);


    %% b) Carga e Potencia Especificada
    Pd   = FC*bus(:, 3)/baseMVA;
    Qd   = FC*bus(:, 4)/baseMVA;
    Sesp = (Pg - Pd) + 1j*(Qg - Qd);


    %% c) Geracao Eolica (opcional)
    if nargin >= 3 && ~isempty(eol)
        for i = 1:numel(eol)
            k = mp(eol(i).barra);

            if tipo(k) ~= 1
                error('Prepara_dados: a barra da fazenda eolica deve ser PQ.');
            end

            Sesp(k) = Sesp(k) + (eol(i).pct/100)*bus(k, 3)/baseMVA;
        end
    end


    %% d) Tap dos Ramos (0 significa ramo sem transformador)
    tap = branch(:, 9);
    tap(tap == 0) = 1;
    branch(:, 9) = tap;



    %% ____________________________________________________________________
    %%
    %%                          3. SAIDA
    %% ____________________________________________________________________
    D.baseMVA = baseMVA;
    D.nb      = nb;
    D.tipo    = tipo;
    D.ref     = ref;
    D.pv      = pv;
    D.pq      = pq;
    D.Vesp    = bus(:, 8);
    D.Sesp    = Sesp;
    D.Qd      = Qd;
    D.Qmin    = Qmin;
    D.Qmax    = Qmax;
    D.bus     = bus;
    D.branch  = branch;
    D.fase0   = branch(:, 10)*pi/180;
    D.num     = num;
    D.kV      = bus(:, 10);

end
