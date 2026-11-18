function [V, a, f, info] = HELMT_Resolver(D, ctrl, lim, metodo, tol, nmax, discreto)
%HELMT_RESOLVER  Laco de eventos comum ao HELM e ao NR com transformadores
%   reguladores: limites de reativo, de tap e de fase, liberacoes e etapa
%   discreta dos taps.
%
%   [V, a, f, info] = HELMT_Resolver(D, ctrl, lim, metodo, tol, nmax, discreto)
%
%   Regras (as mesmas para os dois metodos):
%     - reativo no limite: a barra PV passa a PQ com Q fixo no limite e
%       volta a PV (uma unica vez) se a tensao passar para o lado em que o
%       gerador volta a regular;
%     - tap no limite: o tap fica fixo no limite e e liberado (uma vez) se
%       a tensao controlada passar para o lado regulavel. O sentido e dado
%       pela posicao da barra controlada: no lado do tap (k), aumentar a
%       eleva a tensao; no outro lado (m), reduz;
%     - fase no limite: a fase fica fixa e e liberada (uma vez) se o erro
%       de fluxo passar para o lado regulavel, pelo sinal de dP_km/dphi;
%     - etapa discreta (discreto = true): sem eventos pendentes, cada tap e
%       cada fase livres sao fixados no degrau mais proximo e a solucao e
%       refeita; se a tensao controlada sair da banda morta (lim.banda), o
%       tap anda um degrau no sentido que reduz o erro (no maximo lim.nmov
%       vezes por transformador).
%
%   Entradas
%     D        : estrutura de dados (Prepara_dados)
%     ctrl     : transformadores controlados (ver HELMT_Nucleo)
%     lim      : .amin, .amax, .astep (degrau do tap), .fmin, .fmax,
%                .fstep (degrau da fase, rad), .banda (banda morta de
%                tensao, p.u.), .nmov (movimentos de degrau)
%     metodo   : 'HELM' ou 'NR'
%     tol      : tolerancia (xi do HELM ou desbalanco do NR)
%     nmax     : ordem maxima (HELM) ou iteracoes por passagem (NR)
%     discreto : verdadeiro para a etapa discreta
%
%   Saidas
%     V     : tensoes complexas
%     a, f  : tap e fase finais de cada controle
%     info  : .conv, .passos, .ntot, .red (passagens do HELM aceitas com
%             precisao reduzida, eps < 1e-6), .satQ, .satA, .satF, .disc, .Qger,
%             .registro (texto dos eventos)

    %% ____________________________________________________________________
    %%
    %%                        1. INICIALIZACAO
    %% ____________________________________________________________________

    %% a) Dados
    nb  = D.nb;
    nc  = numel(ctrl);
    pv0 = D.pv;
    tolx = 1e-6;

    if nargin < 7
        discreto = false;
    end


    %% b) Controles: ramos, modos, barra controlada e sentido do tap
    hasV = false(nc, 1);
    hasP = false(nc, 1);
    l1   = zeros(nc, 1);
    sV   = zeros(nc, 1);

    for j = 1:nc
        hasV(j) = any(ctrl(j).modo == 'V');
        hasP(j) = any(ctrl(j).modo == 'P');
        l1(j)   = ctrl(j).ramo(1);

        if hasV(j)
            if ctrl(j).barra == D.branch(l1(j), 1)
                sV(j) = 1;
            else
                sV(j) = -1;
            end
        end
    end


    %% c) Estados dos Eventos
    % satQ: 0 livre, +1 em Qmax, -1 em Qmin; satA/satF: idem para tap/fase
    satQ = zeros(nb, 1);
    satA = zeros(nc, 1);
    satF = zeros(nc, 1);
    libQ = false(nb, 1);
    libA = false(nc, 1);
    libF = false(nc, 1);

    % Valores fixos dos taps e fases (limite ou degrau)
    aFix = D.branch(l1, 9);
    fFix = D.fase0(l1);

    % Etapa discreta: disc = 1 fixado no degrau; movimentos por controle
    disc  = false;
    dA    = false(nc, 1);
    dF    = false(nc, 1);
    nMov  = zeros(nc, 1);

    info.conv     = false;
    info.red      = 0;
    tolred        = 1e-6;
    info.passos   = 0;
    info.ntot     = 0;
    info.registro = {};
    npmax = 30;


    %% ____________________________________________________________________
    %%
    %%                        2. LACO DE EVENTOS
    %% ____________________________________________________________________
    for passo = 1:npmax

        %% a) Caso Modificado pelos Eventos
        D2 = D;
        k  = find(satQ ~= 0);
        D2.tipo(k) = 1;
        Qf = D.Qmin(k);
        Qf(satQ(k) > 0) = D.Qmax(k(satQ(k) > 0));
        D2.Sesp(k) = real(D.Sesp(k)) + 1j*(Qf - D.Qd(k));
        D2.pv = find(D2.tipo == 2);
        D2.pq = find(D2.tipo == 1);

        % i. Controles fixos (limite ou degrau) viram dado; os demais ficam
        cA  = hasV & (satA == 0) & ~dA;
        cF  = hasP & (satF == 0) & ~dF;
        ativ = find(cA | cF);
        ca  = ctrl(ativ);

        for j = 1:nc
            if hasV(j) && ~cA(j)
                D2.branch(ctrl(j).ramo, 9) = aFix(j);
            end

            if hasP(j) && ~cF(j)
                D2.fase0(ctrl(j).ramo) = fFix(j);
            end
        end

        for i = 1:numel(ativ)
            j  = ativ(i);
            md = '';

            if cA(j)
                md = [md, 'V'];
            end

            if cF(j)
                md = [md, 'P'];
            end

            ca(i).modo = md;
        end


        %% b) Solucao
        if strcmpi(metodo, 'HELM')
            [V, at, ft, in] = HELMT_Nucleo(D2, ca, tol, nmax);
            info.ntot = info.ntot + in.n;

            % Passagem aceita com precisao reduzida quando o Wynn estagna
            % abaixo de tolred (as passagens seguintes refazem a solucao)
            if ~in.conv && in.eps < tolred
                in.conv  = true;
                info.red = info.red + 1;
            end
        else
            [V, at, ft, in] = NRT_Nucleo(D2, ca, tol, nmax);
            info.ntot = info.ntot + in.it;
        end

        info.passos = passo;

        if ~in.conv
            break
        end


        %% c) Taps, Fases e Reativo Gerado
        a = aFix;
        f = fFix;

        for i = 1:numel(ativ)
            j = ativ(i);

            if cA(j)
                a(j) = at(i);
            end

            if cF(j)
                f(j) = ft(i);
            end
        end

        lr = arrayfun(@(c) c.ramo(:), ctrl, 'UniformOutput', false);
        lr = vertcat(lr{:}, zeros(0, 1));
        nr = arrayfun(@(c) numel(c.ramo), ctrl);
        ar = repelem_(a, nr);
        fr = repelem_(f, nr);
        Y  = YbusT(D2, lr, ar, fr);
        Qger = imag(V.*conj(Y*V)) + D.Qd;


        %% d) Violacoes de Reativo, Tap e Fase
        ev = false;

        for k = D2.pv(:).'
            if Qger(k) > D.Qmax(k) + tolx
                satQ(k) = 1;
                ev = true;
                info.registro{end+1} = sprintf('Qmax barra %d (passo %d)', k, passo);
            elseif Qger(k) < D.Qmin(k) - tolx
                satQ(k) = -1;
                ev = true;
                info.registro{end+1} = sprintf('Qmin barra %d (passo %d)', k, passo);
            end
        end

        for j = find(cA).'
            if a(j) > lim.amax + 1e-12
                satA(j) = 1;
                aFix(j) = lim.amax;
                ev = true;
            elseif a(j) < lim.amin - 1e-12
                satA(j) = -1;
                aFix(j) = lim.amin;
                ev = true;
            end
        end

        for j = find(cF).'
            if f(j) > lim.fmax + 1e-12
                satF(j) = 1;
                fFix(j) = lim.fmax;
                ev = true;
            elseif f(j) < lim.fmin - 1e-12
                satF(j) = -1;
                fFix(j) = lim.fmin;
                ev = true;
            end
        end

        if ev
            continue
        end


        %% e) Liberacoes (cada unidade no maximo uma vez)
        lib = false;

        for k = find(satQ ~= 0).'
            vk = abs(V(k));

            if ~libQ(k) && ((satQ(k) > 0 && vk > D.Vesp(k) + tolx) || (satQ(k) < 0 && vk < D.Vesp(k) - tolx))
                satQ(k) = 0;
                libQ(k) = true;
                lib     = true;
            end
        end

        for j = find(satA ~= 0).'
            % tap necessario no sentido do erro: da = -sV*(|Vc| - Vesp)
            da = -sV(j)*(abs(V(ctrl(j).barra)) - ctrl(j).Vesp);

            if ~libA(j) && ((satA(j) > 0 && da < -tolx) || (satA(j) < 0 && da > tolx))
                satA(j) = 0;
                libA(j) = true;
                lib     = true;
            end
        end

        for j = find(satF ~= 0).'
            l    = l1(j);
            Pf   = FluxoT(D2, l, V, a(j), f(j));
            sens = sign(FluxoT(D2, l, V, a(j), f(j) + 1e-6) - Pf);
            df   = -(Pf - ctrl(j).Psp)/sens;

            if ~libF(j) && ((satF(j) > 0 && df < -tolx) || (satF(j) < 0 && df > tolx))
                satF(j) = 0;
                libF(j) = true;
                lib     = true;
            end
        end

        if lib
            continue
        end


        %% f) Etapa Discreta
        if discreto && ~disc

            % i. Fixa os taps e fases livres no degrau mais proximo
            for j = find(cA).'
                aFix(j) = Degrau(a(j), aFix(j), lim.astep, lim.amin, lim.amax, D.branch(l1(j), 9));
                dA(j)   = true;
            end

            for j = find(cF).'
                fFix(j) = Degrau(f(j), fFix(j), lim.fstep, lim.fmin, lim.fmax, D.fase0(l1(j)));
                dF(j)   = true;
            end

            disc = true;
            continue
        end

        if disc

            % ii. Banda morta: anda um degrau no sentido que reduz o erro
            mov = false;

            for j = find(dA & satA == 0).'
                err = abs(V(ctrl(j).barra)) - ctrl(j).Vesp;

                if abs(err) > lim.banda && nMov(j) < lim.nmov
                    an = aFix(j) - sV(j)*sign(err)*lim.astep;

                    if an >= lim.amin - 1e-12 && an <= lim.amax + 1e-12
                        aFix(j) = an;
                        nMov(j) = nMov(j) + 1;
                        mov     = true;
                    end
                end
            end

            if mov
                continue
            end
        end

        info.conv = true;
        break
    end


    %% ____________________________________________________________________
    %%
    %%                            3. SAIDA
    %% ____________________________________________________________________
    if ~exist('a', 'var')
        a = aFix;
        f = fFix;
    end

    if ~exist('Qger', 'var')
        Qger = nan(nb, 1);
    end

    info.satQ = satQ;
    info.satA = satA;
    info.satF = satF;
    info.disc = dA | dF;
    info.nmov = nMov;
    info.Qger = Qger;
    info.pv0  = pv0;

end


function x = Degrau(xc, ~, passo, xmin, xmax, x0)
%DEGRAU  Valor discreto mais proximo de xc na grade x0 + i*passo, dentro
%   de [xmin, xmax].
    if passo <= 0
        x = min(max(xc, xmin), xmax);
        return
    end

    x = x0 + round((xc - x0)/passo)*passo;

    while x > xmax + 1e-12
        x = x - passo;
    end

    while x < xmin - 1e-12
        x = x + passo;
    end
end


function y = repelem_(x, n)
%REPELEM_  Repete x(j) n(j) vezes (compativel com controles vazios).
    y = zeros(sum(n), 1);
    fim = cumsum(n(:));

    for j = 1:numel(n)
        y(fim(j)-n(j)+1:fim(j)) = x(j);
    end
end
