function [V, a, fase, info] = NRT_Nucleo(D, ctrl, tol, itmax)
%NRT_NUCLEO  Newton-Raphson completo (esparso) com o modulo do tap e/ou a
%   fase dos transformadores reguladores como incognitas.
%
%   [V, a, fase, info] = NRT_Nucleo(D, ctrl, tol, itmax)
%
%   Incognitas: angulos (exceto V-theta), modulos das barras PQ, modulos de
%   tap (modos V e VP) e fases (modos P e VP). Equacoes: P (exceto V-theta),
%   Q (barras PQ), |Vc| = Vesp (modos V e VP) e P_km = Psp (modos P e VP).
%   Partida plana: |V| = 1 nas barras PQ, |V| = Vesp nas PV e V-theta,
%   angulos nulos, taps e fases do caso.
%
%   Saidas
%     V       : tensoes complexas
%     a, fase : modulo do tap e fase [rad] dos transformadores de ctrl
%     info    : .conv, .it (iteracoes), .erro (maior desbalanco)

    %% ____________________________________________________________________
    %%
    %%                     1. DADOS E INICIALIZACAO
    %% ____________________________________________________________________

    %% a) Dados Gerais
    nb  = D.nb;
    ref = D.ref;
    pv  = D.pv;
    pq  = D.pq;
    nc  = numel(ctrl);
    nth = sort([pv; pq]);
    npq = numel(pq);
    nt  = numel(nth);


    %% b) Controles e Ramos (ctrl(j).ramo pode listar ramos em paralelo no
    %     modo V, com o mesmo tap)
    cb    = zeros(nc, 1);
    Vc    = zeros(nc, 1);
    Psp   = zeros(nc, 1);
    hasV  = false(nc, 1);
    hasP  = false(nc, 1);

    for j = 1:nc
        hasV(j)  = any(ctrl(j).modo == 'V');
        hasP(j)  = any(ctrl(j).modo == 'P');

        if hasV(j)
            cb(j) = ctrl(j).barra;
            Vc(j) = ctrl(j).Vesp;
        end

        if hasP(j)
            Psp(j) = ctrl(j).Psp;
        end
    end

    lst  = arrayfun(@(c) c.ramo(:), ctrl, 'UniformOutput', false);
    nr_  = cellfun(@numel, lst);
    lst  = vertcat(lst{:}, zeros(0, 1));
    dono = zeros(sum(nr_), 1);
    fim  = cumsum(nr_(:));

    for j = 1:nc
        dono(fim(j)-nr_(j)+1:fim(j)) = j;
    end

    u1   = fim - nr_(:) + 1;
    a    = D.branch(lst(u1), 9);
    fase = D.fase0(lst(u1));
    iV   = find(hasV);
    iP   = find(hasP);
    nV   = numel(iV);
    nP   = numel(iP);
    kf   = D.branch(lst, 1);
    km   = D.branch(lst, 2);
    yb   = 1./(D.branch(lst, 3) + 1j*D.branch(lst, 4));
    ykk  = yb + 1j*D.branch(lst, 5)/2;
    pV   = zeros(nc, 1);
    pV(iV) = 1:nV;
    pP   = zeros(nc, 1);
    pP(iP) = 1:nP;
    uP   = u1(iP);


    %% c) Partida Plana e Indices
    Vm      = ones(nb, 1);
    Vm(pv)  = D.Vesp(pv);
    Vm(ref) = D.Vesp(ref);
    Th      = zeros(nb, 1);
    posTh   = zeros(nb, 1);
    posTh(nth) = 1:nt;
    posVm   = zeros(nb, 1);
    posVm(pq)  = nt + (1:npq);
    nx      = nt + npq + nV + nP;
    info.conv = false;
    info.it   = 0;
    info.erro = inf;


    %% ____________________________________________________________________
    %%
    %%                     2. ITERACOES DE NEWTON
    %% ____________________________________________________________________
    for it = 0:itmax

        %% a) Desbalancos
        V  = Vm.*exp(1j*Th);
        ar = a(dono);
        fr = fase(dono);
        Y  = YbusT(D, lst, ar, fr);
        Ib = Y*V;
        S  = V.*conj(Ib);

        % i. Corrente de cada ramo controlado no lado k e fluxo P_km
        T   = ar.*exp(1j*fr);
        Ykk = ykk./ar.^2;
        Ykm = -yb./conj(T);
        Ymk = -yb./T;
        Ik  = Ykk.*V(kf) + Ykm.*V(km);
        Pkm = real(V(kf).*conj(Ik));

        F = [real(S(nth)) - real(D.Sesp(nth)); ...
             imag(S(pq)) - imag(D.Sesp(pq)); ...
             Vm(cb(iV)) - Vc(iV); ...
             Pkm(uP) - Psp(iP)];
        info.erro = max(abs(F));
        info.it   = it;

        if info.erro < tol
            info.conv = true;
            break
        end

        if it == itmax || ~all(isfinite(F))
            break
        end


        %% b) Derivadas das Injecoes em Relacao a V (forma polar, esparsa)
        dV    = spdiags(V./abs(V), 0, nb, nb);
        dgV   = spdiags(V, 0, nb, nb);
        dgI   = spdiags(Ib, 0, nb, nb);
        dSdVm = dgV*conj(Y*dV) + conj(dgI)*dV;
        dSdTh = 1j*dgV*conj(dgI - Y*dgV);


        %% c) Derivadas em Relacao ao Tap e a Fase (so as barras k e m)
        % Ykk = ykk/a^2, Ykm = -y*e^(j*fase)/a, Ymk = -y*e^(-j*fase)/a; os
        % ramos de um grupo somam-se na coluna do seu controle.
        dKa = V(kf).*conj(-2*Ykk./ar.*V(kf) - Ykm./ar.*V(km));
        dMa = V(km).*conj(-Ymk./ar.*V(kf));
        dKf = V(kf).*conj(1j*Ykm.*V(km));
        dMf = V(km).*conj(-1j*Ymk.*V(kf));
        uV  = find(hasV(dono));
        uF  = find(hasP(dono));
        dSa = sparse([kf(uV); km(uV)], [pV(dono(uV)); pV(dono(uV))], [dKa(uV); dMa(uV)], nb, nV);
        dSf = sparse([kf(uF); km(uF)], [pP(dono(uF)); pP(dono(uF))], [dKf(uF); dMf(uF)], nb, nP);


        %% d) Derivadas do Fluxo P_km (ramo unico de cada controle de fluxo)
        % P_km = Re( Vk * conj(Ykk*Vk + Ykm*Vm) )
        k_  = kf(uP);
        m_  = km(uP);
        Vk  = V(k_);
        Vmm = V(m_);
        Ikp = Ik(uP);
        Ykp = Ykk(uP);
        Ymp = Ykm(uP);
        dPth_k = real(1j*Vk.*conj(Ikp) - 1j*Vk.*conj(Ykp).*conj(Vk));
        dPth_m = real(-1j*Vk.*conj(Ymp).*conj(Vmm));
        dPvm_k = real((Vk./abs(Vk)).*conj(Ikp) + Vk.*conj(Ykp.*Vk./abs(Vk)));
        dPvm_m = real(Vk.*conj(Ymp.*Vmm./abs(Vmm)));
        dPf    = real(-1j*Vk.*conj(Ymp).*conj(Vmm));
        dPa    = real(Vk.*conj(-2*Ykp./ar(uP).*Vk - Ymp./ar(uP).*Vmm));


        %% e) Montagem do Jacobiano
        J11 = [real(dSdTh(nth, nth)), real(dSdVm(nth, pq)); ...
               imag(dSdTh(pq, nth)),  imag(dSdVm(pq, pq))];
        J12 = [real(dSa(nth, :)), real(dSf(nth, :)); ...
               imag(dSa(pq, :)),  imag(dSf(pq, :))];

        % i. Linhas de tensao: d|Vc|/d|Vc| = 1
        J21 = sparse(1:nV, posVm(cb(iV)), 1, nV, nt + npq);
        J22 = sparse(nV, nV + nP);

        % ii. Linhas de fluxo (so as variaveis de k e m existentes)
        lin = [];
        col = [];
        val = [];
        bk  = [k_; m_];
        vth = [dPth_k; dPth_m];
        vvm = [dPvm_k; dPvm_m];
        ii  = [(1:nP).'; (1:nP).'];
        ok  = posTh(bk) > 0;
        lin = [lin; ii(ok)];
        col = [col; posTh(bk(ok))];
        val = [val; vth(ok)];
        ok  = posVm(bk) > 0;
        lin = [lin; ii(ok)];
        col = [col; posVm(bk(ok))];
        val = [val; vvm(ok)];
        J31 = sparse(lin, col, val, nP, nt + npq);

        % Fase propria e tap do mesmo controle (modo VP)
        cP  = (1:nP).';
        J32 = sparse(cP, nV + cP, dPf, nP, nV + nP);
        tem = pV(iP) > 0;

        if any(tem)
            J32 = J32 + sparse(cP(tem), pV(iP(tem)), dPa(tem), nP, nV + nP);
        end

        J = [J11, J12; J21, J22; J31, J32];


        %% f) Correcao
        dx = -(J\F);
        Th(nth) = Th(nth) + dx(1:nt);
        Vm(pq)  = Vm(pq) + dx(nt+1:nt+npq);
        a(iV)   = a(iV) + dx(nt+npq+1:nt+npq+nV);
        fase(iP) = fase(iP) + dx(nt+npq+nV+1:nx);
    end

    V = Vm.*exp(1j*Th);

end
