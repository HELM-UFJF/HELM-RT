function [V1, a1, f1, info] = HELMT_Nucleo(D, ctrl, xi, nmax, nestag)
%HELMT_NUCLEO  HELM com transformadores reguladores de relacao complexa
%   como variaveis de estado da recursao (matriz esparsa, Wynn).
%
%   [V1, a1, f1, info] = HELMT_Nucleo(D, ctrl, xi, nmax, nestag)
%
%   Formulacao (mesma do HELM_transformador da tese, com a convencao do
%   MATPOWER para o shunt do lado k, (y + j*b/2)*t^2):
%     - parte nominal do transformador (t_ref = 1/a0, fase de base phi0) na
%       Ybus, separada em parte serie (linhas de soma nula) e shunt;
%     - desvio tau = t - t_ref e phi (fase) como series de Maclaurin, com
%       tau[0] = phi[0] = 0; U = exp(j*phi) pela recursao U' = j*phi'*U;
%     - equacao de controle de tensao |Vc|^2 = 1 + gamma*(Vesp^2 - 1) e de
%       fluxo Re(Vk*conj(Ik)) = gamma*Psp, com Ik igual a corrente que a rede
%       embutida atribui ao ramo.
%   A matriz de recursao e constante, montada por triplas e fatorada uma
%   unica vez. Todas as series (tensoes e estados) sao avaliadas em
%   gamma = 1 pelo algoritmo epsilon de Wynn incremental (Wynn_Passo).
%
%   Entradas
%     D      : estrutura de dados (Prepara_dados)
%     ctrl   : transformadores controlados (.ramo, .modo 'V', 'P' ou 'VP',
%              .barra e .Vesp nos modos V e VP, .Psp nos modos P e VP)
%     xi     : tolerancia de parada (max |estimativa(n) - estimativa(n-1)|)
%     nmax   : ordem maxima
%     nestag : ordens sem melhora ate a parada por estagnacao (padrao 15)
%
%   Saidas
%     V1     : tensoes complexas (p.u.)
%     a1, f1 : modulo do tap e fase (rad) dos transformadores de ctrl
%     info   : .conv, .n (ordem), .eps (ultimo epsilon_n), .singular,
%              .hist (epsilon_n de cada ordem), .histS (maior termo de
%              ordem n, variacao da soma parcial sem acelerar)

    %% a) Argumentos Opcionais
    if nargin < 5 || isempty(nestag)
        nestag = 15;
    end


    %% ____________________________________________________________________
    %%
    %%                1. DADOS E MATRIZ DE RECURSAO
    %% ____________________________________________________________________

    %% a) Dados Gerais
    nb   = D.nb;
    tipo = D.tipo;
    ref  = D.ref;
    pv   = D.pv;
    pq   = D.pq;
    Vesp = D.Vesp;
    Sesp = D.Sesp;
    Pesp = real(Sesp);
    nc   = numel(ctrl);


    %% b) Ybus Nominal: Parte Serie (linhas de soma nula) e Parte Shunt
    Y   = YbusT(D, [], [], []);
    Ysh = full(sum(Y, 2));
    Ysr = Y - spdiags(Ysh, 0, nb, nb);


    %% c) Dados dos Controles e Numeracao dos Estados
    % Os estados de modulo (modo V) vem primeiro e depois os de fase
    % (modo P). A linha de controle de cada estado tem o mesmo indice. No
    % modo V, ctrl(j).ramo pode listar varios ramos em paralelo, que
    % compartilham o mesmo estado tau (mestre-seguidor).
    cb   = zeros(nc, 1);
    Vc   = zeros(nc, 1);
    psp  = zeros(nc, 1);
    hasV = false(nc, 1);
    hasP = false(nc, 1);
    cT   = zeros(nc, 1);
    cP   = zeros(nc, 1);
    uP   = zeros(nc, 1);

    for j = 1:nc
        hasV(j) = any(ctrl(j).modo == 'V');
        hasP(j) = any(ctrl(j).modo == 'P');

        if hasP(j) && numel(ctrl(j).ramo) > 1
            error('HELMT_Nucleo: o controle de fluxo admite um unico ramo.');
        end
    end

    cT(hasV) = 1:nnz(hasV);
    cP(hasP) = nnz(hasV) + (1:nnz(hasP));
    ns       = nnz(hasV) + nnz(hasP);
    iV       = find(hasV);
    iP       = find(hasP);


    %% d) Dados dos Ramos Controlados (um registro por ramo)
    lst  = arrayfun(@(c) c.ramo(:), ctrl, 'UniformOutput', false);
    nr_  = cellfun(@numel, lst);
    lst  = vertcat(lst{:}, zeros(0, 1));
    nu   = numel(lst);
    dono = zeros(sum(nr_), 1);
    fim  = cumsum(nr_(:));

    for j = 1:nc
        dono(fim(j)-nr_(j)+1:fim(j)) = j;
    end

    tref = 1./D.branch(lst, 9);
    yb   = 1./(D.branch(lst, 3) + 1j*D.branch(lst, 4));
    ykk  = yb + 1j*D.branch(lst, 5)/2;
    kf   = D.branch(lst, 1);
    km   = D.branch(lst, 2);
    U0   = exp(1j*D.fase0(lst));
    shf  = ykk.*tref.^2 - yb.*tref.*U0;

    for j = 1:nc
        u = find(dono == j);

        if any(abs(tref(u) - tref(u(1))) > 1e-12)
            error('HELMT_Nucleo: ramos em paralelo devem ter o mesmo tap de base.');
        end

        if hasV(j)
            cb(j) = ctrl(j).barra;
            Vc(j) = ctrl(j).Vesp;

            if tipo(cb(j)) ~= 1
                error('HELMT_Nucleo: a barra controlada deve ser PQ.');
            end
        end

        if hasP(j)
            psp(j) = ctrl(j).Psp;
            uP(j)  = u;
        end
    end


    %% e) Triplas da Parte de Rede
    % Incognitas por barra: PQ e V-theta -> (Re V, Im V); PV -> (Q, Im V).
    % Linhas por barra: partes real e imaginaria da equacao de corrente.
    nx = 2*nb + ns;
    lV = 2*nb + cT(iV);
    lP = 2*nb + cP(iP);
    [kk, mm, g] = find(Ysr);
    ok  = (kk ~= ref) & (mm ~= ref);
    kk  = kk(ok);
    mm  = mm(ok);
    g   = g(ok);
    ePV = (tipo(mm) == 2);

    % i. Colunas de Re V (apenas destinos que nao sao PV)
    r1 = [2*kk(~ePV)-1; 2*kk(~ePV)];
    c1 = [2*mm(~ePV)-1; 2*mm(~ePV)-1];
    v1 = [real(g(~ePV)); imag(g(~ePV))];

    % ii. Colunas de Im V (todos os destinos)
    r2 = [2*kk-1; 2*kk];
    c2 = [2*mm; 2*mm];
    v2 = [-imag(g); real(g)];

    % iii. Barra V-theta: duas linhas identidade
    r3 = [2*ref-1; 2*ref];
    c3 = [2*ref-1; 2*ref];
    v3 = [1; 1];

    % iv. Barras PV: coluna 2k-1 e Q[n], que entra como +jQ[n]
    r4 = 2*pv(:);
    c4 = 2*pv(:) - 1;
    v4 = ones(numel(pv), 1);

    R = [r1; r2; r3; r4];
    C = [c1; c2; c3; c4];
    X = [v1; v2; v3; v4];


    %% f) Triplas das Colunas e Linhas dos Estados (borda)
    % i. Colunas dos estados: contribuicao de cada ramo nas barras k e m
    a1k = 2*tref.*ykk - yb.*U0;
    a1m = -yb.*conj(U0);
    bk  = -1j*yb.*tref.*U0;
    bm  = 1j*yb.*tref.*conj(U0);
    uV  = find(hasV(dono));
    uF  = find(hasP(dono));
    clV = 2*nb + cT(dono(uV));
    clF = 2*nb + cP(dono(uF));
    R   = [R; 2*kf(uV)-1; 2*kf(uV); 2*km(uV)-1; 2*km(uV); ...
              2*kf(uF)-1; 2*kf(uF); 2*km(uF)-1; 2*km(uF)];
    C   = [C; clV; clV; clV; clV; clF; clF; clF; clF];
    X   = [X; real(a1k(uV)); imag(a1k(uV)); real(a1m(uV)); imag(a1m(uV)); ...
              real(bk(uF)); imag(bk(uF)); real(bm(uF)); imag(bm(uF))];

    % ii. Linhas de tensao: Re(Vc[n]) = kappa[n]
    R = [R; 2*nb + cT(hasV)];
    C = [C; 2*cb(hasV) - 1];
    X = [X; ones(nnz(hasV), 1)];

    % iii. Linhas de fluxo (ramo unico)
    for j = find(hasP).'

        % Fluxo: Re(Ik[n]) = ... (ramo unico)
        if hasP(j)
            u    = uP(j);
            col  = 2*nb + cP(j);
            coef = yb(u)*tref(u)*U0(u);
            R    = [R; col];
            C    = [C; col];
            X    = [X; real(-1j*coef)];
            bs   = [kf(u), km(u)];
            sn   = [1, -1];

            for i = 1:2
                b_ = bs(i);

                if b_ == ref
                    continue
                end

                if tipo(b_) ~= 2
                    R = [R; col];
                    C = [C; 2*b_-1];
                    X = [X; sn(i)*real(coef)];
                end

                R = [R; col];
                C = [C; 2*b_];
                X = [X; -sn(i)*imag(coef)];
            end

            % Parcela do tau[n] na linha de fluxo (modo VP)
            if hasV(j)
                R = [R; col];
                C = [C; 2*nb + cT(j)];
                X = [X; real(2*tref(u)*ykk(u) - yb(u)*U0(u))];
            end
        end
    end

    M = sparse(R, C, X, nx, nx);


    %% g) Fatoracao LU (uma unica vez) e Teste de Singularidade
    % Matriz singular: os controles nao sao independentes (por exemplo,
    % fluxo controlado em ramo radial); o problema e mal posto.
    [Lm, Um, Pm, Qm] = lu(M);
    dU = abs(diag(Um));
    info.singular = (min(dU) < 1e-13*max(dU));

    if info.singular
        V1        = nan(nb, 1);
        a1        = nan(nc, 1);
        f1        = nan(nc, 1);
        info.conv = false;
        info.n    = 0;
        info.eps  = inf;
        info.hist  = zeros(0, 1);
        info.histS = zeros(0, 1);
        return
    end


    %% ____________________________________________________________________
    %%
    %%                        2. RECURSAO
    %% ____________________________________________________________________

    %% a) Alocacao e Solucao Germinal (n = 0)
    % Alocacao em blocos de 60 ordens, ampliada quando necessario
    nb_  = min(nmax, 60) + 1;
    V    = zeros(nb, nb_);
    W    = zeros(nb, nb_);
    Q    = zeros(nb, nb_);
    tau  = zeros(nc, nb_);
    phi  = zeros(nc, nb_);
    U    = zeros(nc, nb_);
    Aser = zeros(nu, nb_);
    Bser = zeros(nu, nb_);
    Iks  = zeros(nu, nb_);
    V(:, 1) = 1;
    W(:, 1) = 1;
    U(:, 1) = 1;


    %% b) Wynn: Somas Parciais Empilhadas (tensoes e estados)
    sS       = [V(:, 1); tau(:, 1); phi(:, 1)];
    [wS, ~]  = Wynn_Passo([], sS);
    est_ant  = sS;
    melhor   = sS;
    eps_melh = inf;
    nbest    = 0;
    info.conv = false;
    nfin     = nmax;
    hist     = nan(nmax, 1);
    histS    = nan(nmax, 1);


    %% c) Laco das Ordens n = 1, 2, ...
    for n = 1:nmax

        %% c0) Ampliacao das Matrizes de Coeficientes
        if n + 1 > size(V, 2)
            ac   = min(nmax + 1, 2*size(V, 2)) - size(V, 2);
            V    = [V,    zeros(nb, ac)];
            W    = [W,    zeros(nb, ac)];
            Q    = [Q,    zeros(nb, ac)];
            tau  = [tau,  zeros(nc, ac)];
            phi  = [phi,  zeros(nc, ac)];
            U    = [U,    zeros(nc, ac)];
            Aser = [Aser, zeros(nu, ac)];
            Bser = [Bser, zeros(nu, ac)];
            Iks  = [Iks,  zeros(nu, ac)];
        end


        %% c1) Valores Conhecidos: V-theta (n = 1) e Re V das Barras PV
        K = zeros(nb, 1);

        if n == 1
            K(ref) = Vesp(ref) - 1;
            K(pv)  = (Vesp(pv).^2 - 1)/2;
        else
            K(pv) = -0.5*real(sum(V(pv, 2:n).*conj(V(pv, n:-1:2)), 2));
        end


        %% c2) Lado Direito das Equacoes de Corrente
        r  = zeros(nb, 1);
        Wc = conj(W(:, n));
        r(pq) = conj(Sesp(pq)).*Wc(pq) - Ysh(pq).*V(pq, n);

        if ~isempty(pv)
            Qs = zeros(numel(pv), 1);

            if n > 1
                Qs = sum(Q(pv, 2:n).*conj(W(pv, n:-1:2)), 2);
            end

            r(pv) = Pesp(pv).*Wc(pv) - 1j*Qs - Ysh(pv).*V(pv, n);
        end

        r = r - Ysr*K;


        %% c3) Parte Conhecida das Correntes de Desvio
        % A[n] e B[n] com tau[n] = phi[n] = 0; as parcelas de tau[n] e phi[n]
        % estao na matriz.
        Rn = zeros(nu, 1);
        Un = zeros(nc, 1);

        % i. Parte conhecida de U[n] (por controle)
        if n > 1 && ~isempty(iP)
            Un(iP) = (1j/n)*sum((1:n-1).*phi(iP, 2:n).*U(iP, n:-1:2), 2);
        end

        % ii. Parte conhecida de A[n] e B[n] (por ramo)
        Aser(:, n+1) = sum(tau(dono, 2:n).*tau(dono, n:-1:2), 2);
        Bser(:, n+1) = sum(tau(dono, 2:n).*U(dono, n:-1:2), 2) + tref.*Un(dono);

        % iii. Correntes de desvio conhecidas (convolucoes de ordem n)
        devk = ykk.*sum(Aser(:, 2:n+1).*V(kf, n:-1:1), 2) ...
            - yb.*U0.*sum(Bser(:, 2:n+1).*V(km, n:-1:1), 2);
        devm = -yb.*conj(U0).*sum(conj(Bser(:, 2:n+1)).*V(kf, n:-1:1), 2);
        r    = r - accumarray([kf; km], [devk; devm], [nb, 1]);

        % iv. Parte conhecida de Ik[n]
        Rn = shf.*V(kf, n) + devk;


        %% c4) Vetor do Lado Direito do Sistema Real
        rhs = zeros(nx, 1);
        rhs(1:2:2*nb) = real(r);
        rhs(2:2:2*nb) = imag(r);
        rhs(2*ref-1)  = real(K(ref));
        rhs(2*ref)    = 0;

        % i. Equacoes de tensao: Re(Vc[n]) = kappa[n] (vetorizado)
        if n == 1
            rhs(lV) = (Vc(iV).^2 - 1)/2;
        else
            rhs(lV) = -0.5*real(sum(V(cb(iV), 2:n).*conj(V(cb(iV), n:-1:2)), 2));
        end

        for j = iP.'

            % ii. Equacao de fluxo: Re(Ik[n]) = Psp*delta - sum Re(Vk*conj(Ik))
            if hasP(j)
                u    = uP(j);
                acc  = real(sum(V(kf(u), n:-1:2).*conj(Iks(u, 2:n))));
                alvo = psp(j)*(n == 1);
                lado = alvo - acc - real(Rn(u));
                bs   = [kf(u), km(u)];
                sn   = [1, -1];

                for i = 1:2
                    if tipo(bs(i)) == 2 || bs(i) == ref
                        lado = lado - sn(i)*real(yb(u)*tref(u)*U0(u)*K(bs(i)));
                    end
                end

                rhs(2*nb+cP(j)) = lado;
            end
        end


        %% c5) Solucao (matriz constante ja fatorada)
        x = Qm*(Um\(Lm\(Pm*rhs)));


        %% c6) Coeficientes de Ordem n
        Vn       = x(1:2:2*nb) + 1j*x(2:2:2*nb);
        Vn(pv)   = K(pv) + 1j*x(2*pv);
        Vn(ref)  = K(ref);
        V(:, n+1) = Vn;
        Q(pv, n+1) = x(2*pv-1);

        tau(iV, n+1) = x(lV);
        phi(iP, n+1) = x(lP);
        U(:, n+1)    = Un + 1j*phi(:, n+1);

        Aser(:, n+1) = Aser(:, n+1) + 2*tref.*tau(dono, n+1);
        Bser(:, n+1) = Bser(:, n+1) + tau(dono, n+1) + 1j*tref.*phi(dono, n+1);
        Iks(:, n+1)  = Rn + yb.*tref.*U0.*(V(kf, n+1) - V(km, n+1)) ...
            + (2*tref.*ykk - yb.*U0).*tau(dono, n+1) - 1j*yb.*tref.*U0.*phi(dono, n+1);

        % Recursao do inverso da tensao
        W(:, n+1) = -sum(W(:, 1:n).*V(:, n+1:-1:2), 2);


        %% c7) Wynn e Criterio de Parada
        dS       = [V(:, n+1); tau(:, n+1); phi(:, n+1)];
        sS       = sS + dS;
        histS(n) = max(abs(dS));
        [wS, eS] = Wynn_Passo(wS, sS);
        epsn     = max(abs(eS - est_ant));
        est_ant  = eS;
        hist(n)  = epsn;

        if all(isfinite(eS)) && epsn < eps_melh
            eps_melh = epsn;
            melhor   = eS;
            nbest    = n;
        end

        if epsn < xi
            info.conv = true;
            nfin      = n;
            break
        elseif n - nbest >= nestag
            nfin = n;
            break
        end
    end


    %% ____________________________________________________________________
    %%
    %%                          3. SAIDA
    %% ____________________________________________________________________
    V1  = melhor(1:nb);
    u1  = fim - nr_(:) + 1;
    tp  = tref(u1) + real(melhor(nb+1:nb+nc));
    a1  = 1./tp;
    f1  = angle(U0(u1)) + real(melhor(nb+nc+1:nb+2*nc));

    % Tensao especificada imposta nas barras PV e V-theta
    V1(ref) = Vesp(ref);
    V1(pv)  = Vesp(pv).*V1(pv)./abs(V1(pv));

    info.n    = nfin;
    info.eps  = eps_melh;
    info.hist  = hist(1:nfin);
    info.histS = histS(1:nfin);

end
