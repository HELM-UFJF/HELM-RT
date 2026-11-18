function [P, Q] = FluxoT(D, l, V, a, fase)
%FLUXOT  Fluxo de potencia no lado k (DE) do ramo l, com tap complexo
%   T = a*exp(j*fase) e a convencao do MATPOWER:
%       Ykk = (y + j*b/2)/a^2,  Ykm = -y/conj(T),  S_km = Vk*conj(Ykk*Vk + Ykm*Vm)
%
%   [P, Q] = FluxoT(D, l, V, a, fase)

    %% a) Admitancias do Ramo
    k   = D.branch(l, 1);
    m   = D.branch(l, 2);
    y   = 1./(D.branch(l, 3) + 1j*D.branch(l, 4));
    T   = a.*exp(1j*fase);
    Ykk = (y + 1j*D.branch(l, 5)/2)./a.^2;
    Ykm = -y./conj(T);


    %% b) Fluxo
    S = V(k).*conj(Ykk.*V(k) + Ykm.*V(m));
    P = real(S);
    Q = imag(S);

end
