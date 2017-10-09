function T_init = get_ice_T_init()

global rho g n kc Cp SPY ...
        M M_s N xi dx dzeta dzetadx dzetadx_s zeta hB hB_s hS H H_s dhSdx dhSdx_s dt ...
        de0 Sigma0

LT1 = zeros(N,1);
LT2 = zeros(N,1);
LT3 = zeros(N,1);

T_init = zeros(M,N);

Ts = 268.15 + (5800-hS)*0.007;

Ts(Ts>273.15) = 273.15;
    % be sure the surface ice temperature must not exceed 0 degree
G = 20/1000;
    % geothermal heat flux [W/m2]

for i = 1:M


    RT = zeros(1,N)';

    for j = 2:N-1

        LT3(j) = kc/(rho*Cp*H(i)^2*dzeta^2);
        LT2(j) = -2*kc/(rho*Cp*H(i)^2*dzeta^2);
        LT1(j) = kc/(rho*Cp*H(i)^2*dzeta^2);
    end

    LT1(1) = 0; LT2(1) = -1/dzeta; LT3(1) = 1/dzeta; RT(1) = -G*H(i)/kc;
    LT1(N) = 0; LT2(N) = 1; LT3(N) = 0; RT(N) = Ts(i);

    LT = spdiags([[LT1(2:end);0],LT2,[0;LT3(1:end-1)]],[-1,0,1],N,N);

    T_sol = LT\RT;

%    if T_sol(1) > 273.15 - H(i)*8.7*10^-4
%        LT1(1) = 0; LT2(1) = 1; LT3(1) = 0; 
%        RT(1) = 273.15 - H(i)*8.7*10^-4;
        % prescribed Tpmp bottom BC
%        LT = spdiags([[LT1(2:end);0],LT2,[0;LT3(1:end-1)]],[-1,0,1],N,N);
         % 3-banded sparse matrix
%        T_sol = LT\RT;
         % ice temperature solution vector [K]
%    end

    % do not allow ice temperatures above Tpmp:
    T_pmp = 273.15 - (1-zeta).*H(i)*8.7*10^-4;
    % pressure melting point beneath each node [K]
    index_Tpmp1 = T_sol <= T_pmp; index_Tpmp2 = T_sol > T_pmp;
    % identify where nodes do/do not exceed Tpmp [logic]
    T_sol = T_sol.*index_Tpmp1 + T_pmp.*index_Tpmp2;
        
    T_init(i,:) = T_sol';

end



