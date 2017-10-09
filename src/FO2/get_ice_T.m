function [A, T, W] = get_ice_T(u, u_s, w, w_s, heat_s, T_last)

global rho g n kc Cp SPY ...
        M M_s N xi dx dzeta dzetadx dzetadx_s zeta hB hB_s hS H H_s dhSdx dhSdx_s dt ...
        de0 Sigma0

W = zeros(M,N);
W_vs = zeros(M,N-1);
LT1 = zeros(N,1);
LT2 = zeros(N,1);
LT3 = zeros(N,1);

T = zeros(M,N);
A = zeros(M,N);

Ts = 268.15 + (5708-hS)*0.007;
%Ts = 239 + 8e-8*(xi/1000).^3;
%Ts = 223.15+1.67e-2*(xi/1000);
    % surface annual temperature in 2008-2009; on ~5800 m a.s.l, the annual
    % temperature in 2008-2009 was -5.67 degree.
Ts(Ts>273.15) = 273.15;
    % be sure the surface ice temperature must not exceed 0 degree
G = 20/1000;
    % geothermal heat flux [W/m2]
dt = 1;

for i = 1:M

    u_vs(i,:) = (u(i,1:N-1)+u(i,2:N))/2;
    dzetadx_vs(i,:) = (dzetadx(i,1:N-1)+dzetadx(i,2:N))/2;

    W_vs(i,:) = u_vs(i,:).*dzetadx_vs(i,:) + w_s(i,1:end-1)/H(i);

    W(i,2:end-1) = (W_vs(i,1:end-1)+W_vs(i,2:end))/2;

    if i == 1
        RT = -0*T_last(i,:)'/dt -1*heat_s(i,:)'/rho/Cp/SPY;
        %RT =  -0*heat(i,:)'/rho/Cp/SPY;
    else
        RT = -0*T_last(i,:)'/dt + 1*u_s(i,:)'.*(-T(i-1,:)')/dx/SPY - 1*heat_s(i,:)'/rho/Cp/SPY;
        %RT = u_s(i,:)'.*(-T(i-1,:)')/dx/SPY - heat(i,:)'/rho/Cp/SPY;
    end
    
%    w_vs = [(w(i,1:end-1)+w(i,2:end))/2,0]';
%    if i == 1
%        uu = u(i,:)'*0;
%    else
%        uu = u_s(i,:)';
%        size(uu)
%    end
   

%    LT3 = kc/(rho*Cp*H(i)^2*dzeta^2) - dzetadx(i,:)'.*uu/(2*dzeta)/SPY - w_vs/(2*H(i)*dzeta)/SPY;
%    LT2 = -2*kc/(rho*Cp*H(i)^2*dzeta^2) - uu/dx/SPY; %LT2 = ones(N,1)*LT2;
%    LT1 = kc/(rho*Cp*H(i)^2*dzeta^2) + dzetadx(i,:)'.*uu/(2*dzeta)/SPY + w_vs/(2*H(i)*dzeta)/SPY;

    for j = 2:N-1

        if W(i,j) <= 0
            %fprintf('j: %f\n',j);
%            LT3(j) = kc/(rho*Cp*H(i)^2*dzeta^2) - (u_vs(i,j)*dzetadx_vs(i,j)+w_s(i,j)/H(i))/(dzeta)/SPY;
%            LT2(j) = -2*kc/(rho*Cp*H(i)^2*dzeta^2) - u_s(i,j)/dx/SPY + (u_vs(i,j)*dzetadx_vs(i,j)+w_s(i,j)/H(i))/(dzeta)/SPY;
%            LT1(j) = kc/(rho*Cp*H(i)^2*dzeta^2) + 0*(u_vs(i,j)*dzetadx_vs(i,j)+w_s(i,j)/H(i))/(2*dzeta)/SPY;
%            %LT3(j) = kc/(rho*Cp*H(i)^2*dzeta^2) - W_vs(i,j)/dzeta/SPY;
            LT3(j) = kc/(rho*Cp*H(i)^2*dzeta^2) - 1*W_vs(i,j)/dzeta/SPY;
            LT2(j) = -0*1/dt -2*kc/(rho*Cp*H(i)^2*dzeta^2) - 1*u_s(i,j)/dx/SPY + W_vs(i,j)/dzeta/SPY*1;
            LT1(j) = kc/(rho*Cp*H(i)^2*dzeta^2);
        elseif W(i,j) > 0
            LT3(j) = kc/(rho*Cp*H(i)^2*dzeta^2);
            LT2(j) = -0*1/dt -2*kc/(rho*Cp*H(i)^2*dzeta^2) - u_s(i,j)/dx/SPY*1 - W_vs(i,j-1)/dzeta/SPY*1;
            LT1(j) = kc/(rho*Cp*H(i)^2*dzeta^2) + W_vs(i,j-1)/dzeta/SPY*1;
        end
%
    end
    
%     LT3 = kc/(rho*Cp*H(i)^2*dzeta^2) - dzetadx(i,:).*u_s(i,:)/(2*dzeta)/SPY - w_s(i,:)/(2*H(i)*dzeta)/SPY;
%     LT2 = -2*kc/(rho*Cp*H(i)^2*dzeta^2) - u_s(i,:)/dx/SPY; %LT2 = ones(N,1)*LT2;
%     LT1 = kc/(rho*Cp*H(i)^2*dzeta^2) + dzetadx(i,:).*u_s(i,:)/(2*dzeta)/SPY + w_s(i,:)/(2*H(i)*dzeta)/SPY;
%     LT3 = LT3';LT2 = LT2';LT1 = LT1';


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
%       T_sol = LT\RT;
         % ice temperature solution vector [K]
%    end

    % do not allow ice temperatures above Tpmp:
    T_pmp = 273.15 - (1-zeta).*H(i)*8.7*10^-4;
    % pressure melting point beneath each node [K]
    index_Tpmp1 = T_sol <= T_pmp; index_Tpmp2 = T_sol > T_pmp;
    % identify where nodes do/do not exceed Tpmp [logic]
    T_sol = T_sol.*index_Tpmp1 + T_pmp.*index_Tpmp2;
        
    T(i,:) = T_sol';

    for j=1:N
        if T(i,j)>263.15
            A(i,j) = 1.916e3*exp(-139000./(8.314*T(i,j)));
        end
        if T(i,j)<=263.15
            A(i,j) = 3.985e-13*exp(-60000./(8.314*T(i,j)));
        end
    end

end



