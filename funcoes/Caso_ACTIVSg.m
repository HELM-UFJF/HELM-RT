function [D, c, info] = Caso_ACTIVSg(hora)
%CASO_ACTIVSG  ACTIVSg2000 com os LTCs dos transformadores de rede.
%
%   [D, c, info] = Caso_ACTIVSg(hora)     hora = 5368 (pico) ou 8504 (leve)
%
%   - As barras PV sem faixa de reativo (Qmin = Qmax, geracao renovavel com
%     fator de potencia unitario) passam a PQ com o reativo fixo.
%   - LTCs: transformadores de rede com o lado de baixa tensao (69 a 161 kV)
%     em barra PQ e o lado de alta sem gerador; o tap fica no lado de baixa
%     (lado k) e regula a sua tensao. Transformadores em paralelo na mesma
%     barra formam um grupo com um unico tap (mestre-seguidor).
%   - Tensao especificada: a da barra controlada no caso base da hora
%     (sem controle), arredondada a 0,001 p.u.
%
%   Saidas
%     D    : estrutura de dados (Prepara_dados), com as barras PV ajustadas
%     c    : controles (um por barra de baixa tensao, modo 'V')
%     info : .nPVfixo (barras PV que passaram a PQ), .nRamos (transformadores
%            controlados), .V0 (tensoes do caso base)

    %% a) Caso e Barras PV sem Faixa de Reativo
    SEP = feval(sprintf('case_ACTIVSg2000_h%d', hora));
    D   = Prepara_dados(SEP, 1);
    z   = D.pv(abs(D.Qmax(D.pv) - D.Qmin(D.pv)) < 1e-6);
    D.tipo(z) = 1;
    D.Sesp(z) = real(D.Sesp(z)) + 1j*(D.Qmax(z) - D.Qd(z));
    D.Qmax(z) = inf;
    D.Qmin(z) = -inf;
    D.pv = find(D.tipo == 2);
    D.pq = find(D.tipo == 1);


    %% b) Transformadores de Rede com LTC
    br  = D.branch;
    kv  = D.kV;
    k   = br(:, 1);
    m   = br(:, 2);
    trf = find(SEP.branch(SEP.branch(:, 11) > 0, 9) ~= 0);
    sel = trf(kv(k(trf)) >= 69 & kv(k(trf)) <= 161 & kv(m(trf)) > kv(k(trf)) ...
              & D.tipo(k(trf)) == 1 & D.tipo(m(trf)) ~= 2);
    lvb = unique(k(sel));


    %% c) Caso Base e Tensoes Especificadas
    e0 = struct('ramo', {}, 'modo', {}, 'barra', {}, 'Vesp', {}, 'Psp', {});
    V0 = NRT_Nucleo(D, e0, 1e-10, 30);
    c  = e0;

    for i = 1:numel(lvb)
        c(i).ramo  = sel(k(sel) == lvb(i)).';
        c(i).modo  = 'V';
        c(i).barra = lvb(i);
        c(i).Vesp  = round(abs(V0(lvb(i)))*1000)/1000;
        c(i).Psp   = 0;
    end

    info.nPVfixo = numel(z);
    info.nRamos  = numel(sel);
    info.V0      = V0;

end
