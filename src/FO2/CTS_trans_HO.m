function [T_sol_above, N_CTS]= CTS_trans_HO(i, T, T_sol_CTS, T_pmp, H_CTS, uc, Ts_CTS, Viscc_CTS, de2c_CTS, wc)

%global rho g n kc Cp dx Nx Ny Nz M N I dzeta dxi x b a W hBc hSc dhSdx Hic dzetadx zeta
global rho  kc Cp M dy dzeta dzetadx

dx = dy; I = M;

rho_w = 1000;
    % density of water [kg/m3]
L = 3.35e5; 
    % water latent heat of fusion [J/kg]
mu_tsd = 0;
    % water content threshold. If water content is > this value, water drain away. 
    % Yet this leave behind another problem that we need a water drainage model then.

%Hi_CTS = H_CTS; 
    % keep the original thickness which is used afterwards.

r_CTS = max(find(T_sol_CTS>=T_pmp)); 
    % Find the CTS row number after first temperature calculation. 
    % Then use the CTS as the "bottom" and solve the temperature distribution above this CTS.
    % And find new CTS (bottom) and do new temperature calculation...
    % The loop is shown below

doCTS = 1;
while doCTS == 1

    N_CTS = r_CTS * dzeta; 
    if N_CTS>=1 || H_CTS<1 % if CTS is on th top, break
        T_sol_above = Ts_CTS; % surface temperature
        break;
    end
    
    zeta_CTS = N_CTS:dzeta:1; zeta_CTS = zeta_CTS';
    
    %Kt = 2*A_CTS/(rho*L)*Hi_CTS^6*(rho*g*dhSdx_CTS)^4;
    %mu = -Kt/(Hi_CTS*b_CTS)*quad('(1-x).^4./x',0.0001,N_CTS);
    % calculate water content (Greve, 1997). The horizontal flux is
    % negleced.
    %mu(mu>mu_tsd) = mu_tsd;

    K_CTS = length(zeta_CTS);

    H_CTS = H_CTS*(1-N_CTS);
    
    w = wc(r_CTS+1:end);
    
    LT3 = kc/(rho*Cp*H_CTS^2*dzeta^2) - dzetadx(r_CTS+1:end,i).*uc(r_CTS+1:end)/(2*dzeta) - w/(2*H_CTS*dzeta);
    LT2 = -2*kc/(rho*Cp*H_CTS^2*dzeta^2) - uc(r_CTS+1:end)/dx; %LT2 = ones(K_CTS,1)*LT2;
    LT1 = kc/(rho*Cp*H_CTS^2*dzeta^2) + dzetadx(r_CTS+1:end,i).*uc(r_CTS+1:end)/(2*dzeta) + w/(2*H_CTS*dzeta);
    if i == 1
       RT = uc(r_CTS+1:end).*(-T(r_CTS+1:end,1))/dx - 4*Viscc_CTS(r_CTS+1:end).*de2c_CTS(r_CTS+1:end)/rho/Cp;
    elseif i==I 
       RT = uc(r_CTS+1:end).*(-T(r_CTS+1:end,I-1))/dx - 4*Viscc_CTS(r_CTS+1:end).*de2c_CTS(r_CTS+1:end)/rho/Cp;
    else
       RT = uc(r_CTS+1:end).*(-T(r_CTS+1:end,i-1))/dx - 4*Viscc_CTS(r_CTS+1:end).*de2c_CTS(r_CTS+1:end)/rho/Cp;
    end
    
    % boundary condition
    LT2(1) = 1; LT3(1) = 0; RT(1) = 273.15-8.7e-4*H_CTS;
    LT1(K_CTS) = 0; LT2(K_CTS) = 1; RT(K_CTS) = Ts_CTS;
    
    LT = spdiags([[LT1(2:end);0],LT2,[0;LT3(1:end-1)]],[-1,0,1],K_CTS,K_CTS);
    
    T_sol_above = LT\RT;
    
    if T_sol_above(2) < T_sol_above(1)+8.7e-4*H_CTS/K_CTS-mu_tsd*rho_w*L*w(1)*H_CTS/kc/K_CTS;
        % Check if the result satisfy the gradient condition. If not, move up 1 node.
        doCTS = 0;
    else
        r_CTS = r_CTS+1;
    end
end