function [w, est] = Wynn_Passo(w, s)
%WYNN_PASSO  Algoritmo epsilon de Wynn incremental (uma soma parcial por vez).
%
%   A cada nova soma parcial s = S_n (uma por linha) calcula-se apenas a
%   nova antidiagonal da tabela epsilon pela regra do losango
%       eps_{j+1}^(k) = eps_{j-1}^(k+1) + 1/(eps_j^(k+1) - eps_j^(k)),
%   em que eps_j^(k+1) esta na antidiagonal nova e eps_j^(k), eps_{j-1}^(k+1)
%   na anterior. O custo por ordem e O(n) por serie, contra O(n^2) da
%   reconstrucao completa da tabela (Wynn_Epsilon), com o mesmo resultado.
%
%   w   : estado ([] na primeira chamada, com s = S_0)
%   est : estimativa do limite = elemento de ordem par mais alta da
%         antidiagonal nova (coluna par mais alta da tabela)


    %% a) Primeira Soma Parcial
    % O estado guarda tambem o modulo dos elementos da antidiagonal (w.mod),
    % calculado uma unica vez, para o teste de quebra da tabela.
    s  = s(:);
    nl = numel(s);

    if isempty(w)
        % Antidiagonal: niveis 0, 1, ...
        w.ant = s;
        w.mod = abs(s);

        % Nivel mais alto valido por linha
        w.dep = zeros(nl, 1);
        est   = s;
        return
    end


    %% b) Nova Antidiagonal
    nAnt = size(w.ant, 2);
    nov  = zeros(nl, nAnt + 1);
    mdN  = zeros(nl, nAnt + 1);
    nov(:,1) = s;
    mdN(:,1) = abs(s);
    dep  = zeros(nl, 1);
    ativa = true(nl, 1);

    for j = 0:nAnt-1
        % Eps_j da antidiagonal anterior existe
        ativa = ativa & (w.dep >= j);

        if ~any(ativa)
            break
        end

        % i. Diferenca e teste de quebra (serie estacionaria)
        d  = nov(:,j+1) - w.ant(:,j+1);
        ad = abs(d);
        ativa = ativa & (ad > 1e-13*max(max(mdN(:,j+1), w.mod(:,j+1)), 1e-300));

        % Valor inocuo
        d(~ativa) = 1;

        % ii. Regra do losango
        if j == 0
            e = 1./d;
        else
            e = w.ant(:,j) + 1./d;
        end

        me    = abs(e);
        ativa = ativa & (me < Inf);
        nov(:,j+2) = e;
        mdN(:,j+2) = me;
        dep(ativa) = j + 1;
    end


    %% c) Estimativa (nivel par mais alto valido)
    par = dep - mod(dep, 2);
    est = nov((1:nl).' + nl*par);

    w.ant = nov;
    w.mod = mdN;
    w.dep = dep;
end
