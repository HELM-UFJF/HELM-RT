function Y = YbusT(D, ramos, a, fase)
%YBUST  Matriz de admitancias nodais esparsa com tap complexo.
%
%   Y = YbusT(D, ramos, a, fase)
%
%   Convencao do MATPOWER: ramo de k (DE) para m (PARA), relacao complexa
%   T = a*exp(j*fase) no lado k:
%       Ykk = (y + j*b/2)/a^2     Ykm = -y/conj(T)
%       Ymk = -y/T                Ymm = y + j*b/2
%
%   Entradas
%     D     : estrutura de dados (Prepara_dados)
%     ramos : ramos com tap e/ou fase alterados ([] = nenhum)
%     a     : novos modulos de tap desses ramos ([] = mantem D.branch)
%     fase  : novas fases [rad] desses ramos ([] = mantem D.fase0)

    %% a) Parametros dos Ramos
    br  = D.branch;
    nb  = D.nb;
    k   = br(:, 1);
    m   = br(:, 2);
    tap = br(:, 9);
    phi = D.fase0(:);

    if ~isempty(ramos)
        if ~isempty(a)
            tap(ramos) = a(:);
        end

        if ~isempty(fase)
            phi(ramos) = fase(:);
        end
    end

    T   = tap.*exp(1j*phi);
    y   = 1./(br(:, 3) + 1j*br(:, 4));
    bsh = 1j*br(:, 5)/2;


    %% b) Elementos dos Ramos
    Ykk = (y + bsh)./(tap.^2);
    Ymm = y + bsh;
    Ykm = -y./conj(T);
    Ymk = -y./T;


    %% c) Montagem Esparsa (com os bancos shunt das barras)
    lin = [k; m; k; m; (1:nb).'];
    col = [k; m; m; k; (1:nb).'];
    val = [Ykk; Ymm; Ykm; Ymk; 1j*D.bus(:, 6)/D.baseMVA];
    Y   = sparse(lin, col, val, nb, nb);

end
