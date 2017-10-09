function [u_s, u] = get_ice_u_periodic(A_s, visc_s, T_s, u_lst);

global rho g n kc Cp SPY ...
        M M_s N xi dx dzeta dzetadx dzetadx_s zeta hB hB_s hS H H_s dhSdx dhSdx_s dt ...
        de0 Sigma0

[LHS, RHS] = build_V_lhs_rhs_periodic(A_s, visc_s, T_s, u_lst);

v = LHS\RHS;
%v = linsolve(LHS,RHS);

for k=1:M_s*N
    if fix(k/N) == k/N
        i = k/N;
    else
        i = fix(k/N) + 1;
    end
    j = k-(i-1)*N;
    u_s(i,j) = v(k);
end

u = [u_s(1,:);(u_s(2:M_s-2,:)+u_s(3:M_s-1,:))/2;u_s(M_s,:)];
