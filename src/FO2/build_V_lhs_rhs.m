function [LHS, RHS] = build_V_lhs_rhs(A_s, visc_s, T_s)

global rho g n kc Cp SPY ...
        M M_s N xi dx dzeta dzetadx dzetadx_s zeta hB hB_s hS H H_s dhSdx dhSdx_s dt ...
        de0 Sigma0

LHS = zeros(M_s*N, M_s*N);
RHS = zeros(M_s*N,1);

%dzetadx_s = zeros(M_s,N);

L1 = zeros(M_s,N);
L2 = zeros(M_s,N);
L3 = zeros(M_s,N);
L4 = zeros(M_s,N);
L5 = zeros(M_s,N);
Ls = zeros(M_s,N);

PHTG = 1;

%--------------------- build Right Hand Side vector --------------------
for i = 1:M_s
for j = 1:N
      
    k = (i-1)*N + j;
    f = [0, 0.313, 0.558, 0.790, 0.884, 0.929, 0.954, 0.990, 1]; % rec
    f1 = [0, 0.251, 0.448, 0.653, 0.748, 0.803, 0.839, 0.917, 1]; % U
    ar = [0, 0.5, 1, 2, 3, 4, 5, 10, 1e4];
    fi = interp1(ar, f, 600/H_s(i));
    RHS(k,1) = 1*rho*g*dhSdx_s(i);
        
end
end
%-------------------------- end RHS -------------------------------------

%------------------- build Left Hand Side matrix ------------------------    
for i = 2:M_s-1
for j = 2:N-1

    k = (i-1)*N + j; % mapping method
    
    dzetadx_s(1,j) = -((hB_s(1+1)-hB_s(1))/(dx/2)+...
                    (j-1)*dzeta*(H_s(1+1)-H_s(1))/(dx/2))/H_s(1);
    if i == 2 && PHTG
    dzetadx_s(i,j) = -1/H_s(i)*(1/(3*dx)*(-4*hB_s(i-1)+3*hB_s(i)+hB_s(i+1)) +...
                    (j-1)*dzeta/(3*dx)*(-4*H_s(i-1)+3*H_s(i)+H_s(i+1)));
    dzetadx_s(i,N) = -1/H_s(i)*(1/(3*dx)*(-4*hB_s(i-1)+3*hB_s(i)+hB_s(i+1)) +...
                    (N-1)*dzeta/(3*dx)*(-4*H_s(i-1)+3*H_s(i)+H_s(i+1)));
    dzetadx_s(i,1) = -1/H_s(i)*(1/(3*dx)*(-4*hB_s(i-1)+3*hB_s(i)+hB_s(i+1)) +...
                    (1-1)*dzeta/(3*dx)*(-4*H_s(i-1)+3*H_s(i)+H_s(i+1)));
    elseif i == M_s-1 && PHTG
    dzetadx_s(i,j) = -1/H_s(i)*(1/(3*dx)*(-hB_s(i-1)-3*hB_s(i)+4*hB_s(i+1)) +...
                    (j-1)*dzeta/(3*dx)*(-H_s(i-1)-3*H_s(i)+4*H_s(i+1)));
    dzetadx_s(i,N) = -1/H_s(i)*(1/(3*dx)*(-hB_s(i-1)-3*hB_s(i)+4*hB_s(i+1)) +...
                    (N-1)*dzeta/(3*dx)*(-H_s(i-1)-3*H_s(i)+4*H_s(i+1)));
    dzetadx_s(i,1) = -1/H_s(i)*(1/(3*dx)*(-hB_s(i-1)-3*hB_s(i)+4*hB_s(i+1)) +...
                    (1-1)*dzeta/(3*dx)*(-H_s(i-1)-3*H_s(i)+4*H_s(i+1)));
    else
    dzetadx_s(i,j) = -((hB_s(i+1)-hB_s(i-1))/(2*dx)+...
                    (j-1)*dzeta*(H_s(i+1)-H_s(i-1))/(2*dx))/H_s(i);
    dzetadx_s(i,N) = -((hB_s(i+1)-hB_s(i-1))/(2*dx)+...
                    (N-1)*dzeta*(H_s(i+1)-H_s(i-1))/(2*dx))/H_s(i);
    dzetadx_s(i,1) = -((hB_s(i+1)-hB_s(i-1))/(2*dx)+...
                    (1-1)*dzeta*(H_s(i+1)-H_s(i-1))/(2*dx))/H_s(i);
    end
    dzetadx_s(M_s,j) = -((hB_s(M_s)-hB_s(M_s-1))/(dx/2)+...
                    (j-1)*dzeta*(H_s(M_s)-H_s(M_s-1))/(dx/2))/H_s(M_s);
    dzetadx_s(M_s,N) = -((hB_s(M_s)-hB_s(M_s-1))/(dx/2)+...
                    (N-1)*dzeta*(H_s(M_s)-H_s(M_s-1))/(dx/2))/H_s(M_s);

    L1(i,j) = 4*visc_s(i,j); 
    L1(i,N) = 4*visc_s(i,N);
    L1(i,1) = 4*visc_s(i,1);
    L1(1,j) = 4*visc_s(1,j);
    

    L2(i,j) = 4*visc_s(i,j)*((dzetadx_s(i,j))^2 + 1/4*(1/H_s(i))^2);
    L2(i,N) = 4*visc_s(i,N)*((dzetadx_s(i,N))^2 + 1/4*(1/H_s(i))^2);
    L2(i,1) = 4*visc_s(i,1)*((dzetadx_s(i,1))^2 + 1/4*(1/H_s(i))^2);
    L2(1,j) = 4*visc_s(1,j)*((dzetadx_s(1,j))^2 + 1/4*(1/H_s(i))^2);

    L3(i,j) = 8*visc_s(i,j)*dzetadx_s(i,j);
    L3(i,N) = 8*visc_s(i,N)*dzetadx_s(i,N);
    L3(i,1) = 8*visc_s(i,1)*dzetadx_s(i,1);
    L3(1,j) = 8*visc_s(1,j)*dzetadx_s(1,j);

     if i == 2 && PHTG
     L4(i,j) = 4*visc_s(i,j)*(1/(3*dx)*(-4*dzetadx_s(i-1,j)+3*dzetadx_s(i,j)+dzetadx_s(i+1,j)))+...
             4*visc_s(i,j)*dzetadx_s(i,j)*(dzetadx_s(i,j+1)-dzetadx_s(i,j-1))/(2*dzeta) +...
             4*dzetadx_s(i,j)*(1/(3*dx)*(-4*visc_s(i-1,j)+3*visc_s(i,j)+visc_s(i+1,j))) +...
             (visc_s(i,j+1)-visc_s(i,j-1))/(2*dzeta)*(4*dzetadx_s(i,j)^2+(1/H_s(i))^2);
     L4(i,N) = 4*visc_s(i,N)*(1/(3*dx)*(-4*dzetadx_s(i-1,N)+3*dzetadx_s(i,N)+dzetadx_s(i+1,N)))+...
             4*visc_s(i,N)*dzetadx_s(i,N)*(dzetadx_s(i,N)-dzetadx_s(i,N-1))/(dzeta) +...
             4*dzetadx_s(i,N)*(1/(3*dx)*(-4*visc_s(i-1,N)+3*visc_s(i,N)+visc_s(i+1,N))) +...
             (visc_s(i,N)-visc_s(i,N-1))/(dzeta)*(4*dzetadx_s(i,N)^2+(1/H_s(i))^2);
    L4(i,1) = 4*visc_s(i,1)*(1/(3*dx)*(-4*dzetadx_s(i-1,1)+3*dzetadx_s(i,1)+dzetadx_s(i+1,1)))+...
            4*visc_s(i,1)*dzetadx_s(i,1)*(dzetadx_s(i,2)-dzetadx_s(i,1))/(dzeta) +...
            4*dzetadx_s(i,1)*(1/(3*dx)*(-4*visc_s(i-1,1)+3*visc_s(i,1)+visc_s(i+1,1))) +...
            (visc_s(i,2)-visc_s(i,1))/(dzeta)*(4*dzetadx_s(i,1)^2+(1/H_s(i))^2);
    elseif i == M_s-1 && PHTG && 1
    L4(i,j) = 4*visc_s(i,j)*(dzetadx_s(i+1,j)-dzetadx_s(i-1,j))/...
                (2*dx)+ 4*visc_s(i,j)*dzetadx_s(i,j)*...
                (dzetadx_s(i,j+1)-dzetadx_s(i,j-1))/(2*dzeta) +...
              4*(visc_s(i+1,j)-visc_s(i-1,j))/(2*dx)*dzetadx_s(i,j) +...
               (visc_s(i,j+1)-visc_s(i,j-1))/(2*dzeta)*(4*dzetadx_s(i,j)^2+(1/H_s(i))^2);
    %L4(i,j) = 3*L4(i,j);
    L41(i,j) = 4*visc_s(i,j)*(1/(3*dx)*(-dzetadx_s(i-1,j)-3*dzetadx_s(i,j)+4*dzetadx_s(i+1,j)))+...
            4*visc_s(i,j)*dzetadx_s(i,j)*(dzetadx_s(i,j+1)-dzetadx_s(i,j-1))/(2*dzeta) +...
            4*dzetadx_s(i,j)*(1/(3*dx)*(-visc_s(i-1,j)-3*visc_s(i,j)+4*visc_s(i+1,j))) +...
            (visc_s(i,j+1)-visc_s(i,j-1))/(2*dzeta)*(4*dzetadx_s(i,j)^2+(1/H_s(i))^2);
    %fprintf('%f, %f, %f\n',L4(i,j),L41(i,j),L41(i,j)/L4(i,j));
    L4(i,N) = 4*visc_s(i,N)*(dzetadx_s(i+1,N)-dzetadx_s(i-1,N))/...
                (2*dx)+ 4*visc_s(i,N)*dzetadx_s(i,N)*...
                (dzetadx_s(i,N)-dzetadx_s(i,N-1))/(dzeta) + ...
              4*(visc_s(i+1,N)-visc_s(i-1,N))/(2*dx)*dzetadx_s(i,N) +...
               (visc_s(i,N)-visc_s(i,N-1))/(dzeta)*(4*dzetadx_s(i,N)^2+(1/H_s(i))^2);
    L41(i,N) = 4*visc_s(i,N)*(1/(3*dx)*(-dzetadx_s(i-1,N)-3*dzetadx_s(i,N)+4*dzetadx_s(i+1,N)))+...
            4*visc_s(i,N)*dzetadx_s(i,N)*(dzetadx_s(i,N)-dzetadx_s(i,N-1))/(dzeta) +...
            4*dzetadx_s(i,N)*(1/(3*dx)*(-visc_s(i-1,N)-3*visc_s(i,N)+4*visc_s(i+1,N))) +...
            (visc_s(i,N)-visc_s(i,N-1))/(dzeta)*(4*dzetadx_s(i,N)^2+(1/H_s(i))^2);
    L4(i,1) = 4*visc_s(i,1)*(1/(3*dx)*(-dzetadx_s(i-1,1)-3*dzetadx_s(i,1)+4*dzetadx_s(i+1,1)))+...
            4*visc_s(i,1)*dzetadx_s(i,1)*(dzetadx_s(i,1+1)-dzetadx_s(i,1))/(dzeta) +...
            4*dzetadx_s(i,1)*(1/(3*dx)*(-visc_s(i-1,1)-3*visc_s(i,1)+4*visc_s(i+1,1))) +...
            (visc_s(i,1+1)-visc_s(i,1))/(dzeta)*(4*dzetadx_s(i,1)^2+(1/H_s(i))^2);
    else
    L4(i,j) = 4*visc_s(i,j)*(dzetadx_s(i+1,j)-dzetadx_s(i-1,j))/...
                (2*dx)+ 4*visc_s(i,j)*dzetadx_s(i,j)*...
                (dzetadx_s(i,j+1)-dzetadx_s(i,j-1))/(2*dzeta) +...
              4*(visc_s(i+1,j)-visc_s(i-1,j))/(2*dx)*dzetadx_s(i,j) +...
               (visc_s(i,j+1)-visc_s(i,j-1))/(2*dzeta)*(4*dzetadx_s(i,j)^2+(1/H_s(i))^2);
    L4(i,N) = 4*visc_s(i,N)*(dzetadx_s(i+1,N)-dzetadx_s(i-1,N))/...
                (2*dx)+ 4*visc_s(i,N)*dzetadx_s(i,N)*...
                (dzetadx_s(i,N)-dzetadx_s(i,N-1))/(dzeta) + ...
              4*(visc_s(i+1,N)-visc_s(i-1,N))/(2*dx)*dzetadx_s(i,N) +...
               (visc_s(i,N)-visc_s(i,N-1))/(dzeta)*(4*dzetadx_s(i,N)^2+(1/H_s(i))^2);
    L4(i,1) = 4*visc_s(i,1)*(dzetadx_s(i+1,1)-dzetadx_s(i-1,1))/(2*dx)+...
                 4*visc_s(i,1)*dzetadx_s(i,1)*(dzetadx_s(i,2)-dzetadx_s(i,1))/(dzeta) + ...
              4*(visc_s(i+1,1)-visc_s(i-1,1))/(2*dx)*dzetadx_s(i,1) +...
               (visc_s(i,2)-visc_s(i,1))/(dzeta)*(4*dzetadx_s(i,1)^2+(1/H_s(i))^2);        
    L4(1,j) = 4*visc_s(1,j)*(dzetadx_s(2,j)-dzetadx_s(1,j))/...
                (dx)+ 4*visc_s(1,j)*dzetadx_s(1,j)*...
                (dzetadx_s(1,j+1)-dzetadx_s(1,j-1))/(2*dzeta) +...
              4*(visc_s(2,j)-visc_s(1,j))/(dx)*dzetadx_s(1,j) +...
               (visc_s(1,j+1)-visc_s(1,j-1))/(2*dzeta)*(4*dzetadx_s(1,j)^2+(1/H_s(1))^2);
    end

     if i == 2 && PHTG
     L5(i,j) = 4*(1/(3*dx)*(-4*visc_s(i-1,j)+3*visc_s(i,j)+visc_s(i+1,j))) +...
             4*(visc_s(i,j+1)-visc_s(i,j-1))/(2*dzeta)*(dzetadx_s(i,j));
     L5(i,N) = 4*(1/(3*dx)*(-4*visc_s(i-1,N)+3*visc_s(i,N)+visc_s(i+1,N))) +...
             4*(visc_s(i,N)-visc_s(i,N-1))/(dzeta)*(dzetadx_s(i,N));
    L5(i,1) = 4*(1/(3*dx)*(-4*visc_s(i-1,1)+3*visc_s(i,j)+visc_s(i+1,1))) +...
            4*(visc_s(i,2)-visc_s(i,1))/(dzeta)*(dzetadx_s(i,1));
    elseif i == M_s-1 && PHTG && 1
    L5(i,j) = 4*(1/(3*dx)*(-visc_s(i-1,j)-3*visc_s(i,j)+4*visc_s(i+1,j))) +...
            4*(visc_s(i,j+1)-visc_s(i,j-1))/(2*dzeta)*(dzetadx_s(i,j));
    L5(i,N) = 4*(1/(3*dx)*(-visc_s(i-1,N)-3*visc_s(i,N)+4*visc_s(i+1,N))) +...
            4*(visc_s(i,N)-visc_s(i,N-1))/(dzeta)*(dzetadx_s(i,N));
    L5(i,1) = 4*(1/(3*dx)*(-visc_s(i-1,1)-3*visc_s(i,1)+4*visc_s(i+1,1))) +...
            4*(visc_s(i,1+1)-visc_s(i,1))/(dzeta)*(dzetadx_s(i,1));
    else
    L5(i,j) = 4*(visc_s(i+1,j)-visc_s(i-1,j))/(2*dx) + ...
        4*(visc_s(i,j+1)-visc_s(i,j-1))/(2*dzeta)*(dzetadx_s(i,j));
    L5(i,N) = 4*(visc_s(i+1,N)-visc_s(i-1,N))/(2*dx) + ...
        4*(visc_s(i,N)-visc_s(i,N-1))/(dzeta)*(dzetadx_s(i,N));
    L5(i,1) = 4*(visc_s(i+1,1)-visc_s(i-1,1))/(2*dx) + ...
        4*(visc_s(i,2)-visc_s(i,1))/(dzeta)*(dzetadx_s(i,1));
    end

    Ls(i,N) = 4*dhSdx_s(i)*dzeta/dx/(1/H_s(i)-4*dhSdx_s(i)*dzetadx_s(i,N));
    Ls(1,N) = 4*dhSdx_s(1)*dzeta/dx/(1/H_s(1)-4*dhSdx_s(1)*dzetadx_s(1,N));
    Ls(M_s,N) = 4*dhSdx_s(M_s)*dzeta/dx/(1/H_s(M_s)-4*dhSdx_s(M_s)*dzetadx_s(M_s,N));

     if i == 2 && PHTG
         LHS(k,k-N) = 8*L1(i,j)/(3*dx^2) - 4*L5(i,j)/(3*dx);
         LHS(k,k) = -4*L1(i,j)/dx^2 - 2*L2(i,j)/dzeta^2 + L5(i,j)/dx;
         LHS(k,k+N) = 4*L1(i,j)/(3*dx^2) + L5(i,j)/(3*dx);
         LHS(k,k+1) = L2(i,j)/dzeta^2 + L3(i,j)/(2*dx*dzeta) + L4(i,j)/(2*dzeta);
         LHS(k,k-1) = L2(i,j)/dzeta^2 - L3(i,j)/(2*dx*dzeta) - L4(i,j)/(2*dzeta);
         LHS(k,k-N+1) = -2*L3(i,j)/(3*dx*dzeta);
         LHS(k,k-N-1) = 2*L3(i,j)/(3*dx*dzeta);
         LHS(k,k+N+1) = L3(i,j)/(6*dx*dzeta);
         LHS(k,k+N-1) = -L3(i,j)/(6*dx*dzeta);
    elseif i == M_s-1 && PHTG && 1
        LHS(k,k-N) = 4*L1(i,j)/(3*dx^2) - L5(i,j)/(3*dx);
        LHS(k,k) = -4*L1(i,j)/dx^2 - 2*L2(i,j)/dzeta^2 - L5(i,j)/dx;
        LHS(k,k+N) = 8*L1(i,j)/(3*dx^2) + 4*L5(i,j)/(3*dx);
        LHS(k,k+1) = L2(i,j)/dzeta^2 - L3(i,j)/(2*dx*dzeta) + L4(i,j)/(2*dzeta);
        LHS(k,k-1) = L2(i,j)/dzeta^2 + L3(i,j)/(2*dx*dzeta) - L4(i,j)/(2*dzeta);
        LHS(k,k-N+1) = -L3(i,j)/(6*dx*dzeta);
        LHS(k,k-N-1) = L3(i,j)/(6*dx*dzeta);
        LHS(k,k+N+1) = 2*L3(i,j)/(3*dx*dzeta);
        LHS(k,k+N-1) = -2*L3(i,j)/(3*dx*dzeta);
    else
        LHS(k,k-N-1) = L3(i,j)/(4*dx*dzeta);
        LHS(k,k-N) = L1(i,j)/dx^2 - L5(i,j)/(2*dx);
        LHS(k,k-N+1) = -L3(i,j)/(4*dx*dzeta);
        LHS(k,k-1) = L2(i,j)/dzeta^2 - L4(i,j)/(2*dzeta);
        LHS(k,k) = -2*L1(i,j)/dx^2 - 2*L2(i,j)/dzeta^2;
        LHS(k,k+1) = L2(i,j)/dzeta^2 + L4(i,j)/(2*dzeta);
        LHS(k,k+N-1) = -L3(i,j)/(4*dx*dzeta);
        LHS(k,k+N) = L1(i,j)/dx^2 + L5(i,j)/(2*dx);
        LHS(k,k+N+1) = L3(i,j)/(4*dx*dzeta);
    end

end
end
%------------------------------ end LHS ------------------------------------

%-------------------------- implement boundary condition -------------------
% 1) surface: when i=2:M_s-1, j=N
for i = 2:M_s-1
 if i == 2 && PHTG
     LHS(i*N,(i-1)*N) = 8*L1(i,N)/(3*dx^2) - 4*L5(i,N)/(3*dx) -...
                     8*Ls(i,N)/3*(L2(i,N)/dzeta^2 + L3(i,N)/(2*dx*dzeta) + L4(i,N)/(2*dzeta)) +...
                     4*L3(i,N)/(3*dx*dzeta)*Ls(i-1,N);
     LHS(i*N,i*N) = -4*L1(i,N)/dx^2 - 2*L2(i,N)/dzeta^2 + L5(i,N)/dx +...
                 2*Ls(i,N)*(L2(i,N)/dzeta^2 + L3(i,N)/(2*dx*dzeta) + L4(i,N)/(2*dzeta)) - ...
                 4*L3(i,N)/(9*dx*dzeta)*Ls(i-1,N) - 4*L3(i,N)/(9*dx*dzeta)*Ls(i+1,N);
     LHS(i*N,(i+1)*N) = 4*L1(i,N)/(3*dx^2) + L5(i,N)/(3*dx) +...
                     2*Ls(i,N)/3*(L2(i,N)/dzeta^2 + L3(i,N)/(2*dx*dzeta) + L4(i,N)/(2*dzeta)) +...
                     L3(i,N)/(3*dx*dzeta)*Ls(i+1,N);
     LHS(i*N,i*N-1) = L2(i,N)/dzeta^2 - L3(i,N)/(2*dx*dzeta) - L4(i,N)/(2*dzeta) +...
                     (L2(i,N)/dzeta^2 + L3(i,N)/(2*dx*dzeta) + L4(i,N)/(2*dzeta));
     LHS(i*N,(i+2)*N) = L3(i,N)/(9*dx*dzeta)*Ls(i+1,N);
elseif i == M_s-1 && PHTG && 1
    LHS(i*N,(i-1)*N) = 4*L1(i,N)/(3*dx^2) - L5(i,N)/(3*dx) -...
                    2*Ls(i,N)/3*(L2(i,N)/dzeta^2 - L3(i,N)/(2*dx*dzeta) + L4(i,N)/(2*dzeta)) +...
                    L3(i,N)/(3*dx*dzeta)*Ls(i-1,N);
    LHS(i*N,i*N) = -4*L1(i,N)/dx^2 - 2*L2(i,N)/dzeta^2 - L5(i,N)/dx -...
                2*Ls(i,N)*(L2(i,N)/dzeta^2 - L3(i,N)/(2*dx*dzeta) + L4(i,N)/(2*dzeta)) - ...
                4*L3(i,N)/(9*dx*dzeta)*Ls(i-1,N) - 4*L3(i,N)/(9*dx*dzeta)*Ls(i+1,N);
    LHS(i*N,(i+1)*N) = 8*L1(i,N)/(3*dx^2) + 4*L5(i,N)/(3*dx) +...
                    8*Ls(i,N)/3*(L2(i,N)/dzeta^2 - L3(i,N)/(2*dx*dzeta) + L4(i,N)/(2*dzeta)) -...
                    4*L3(i,N)/(3*dx*dzeta)*Ls(i+1,N);
    LHS(i*N,i*N-1) = L2(i,N)/dzeta^2 + L3(i,N)/(2*dx*dzeta) - L4(i,N)/(2*dzeta) +...
                    L2(i,N)/dzeta^2 - L3(i,N)/(2*dx*dzeta) + L4(i,N)/(2*dzeta);
    LHS(i*N,(i-2)*N) = L3(i,N)/(9*dx*dzeta)*Ls(i-1,N);
else
    LHS(i*N,(i+1)*N) = L1(i,N)/dx^2 + L2(i,N)*Ls(i,N)/dzeta^2 +...
                     L4(i,N)*Ls(i,N)/(2*dzeta) + L5(i,N)/(2*dx);
    LHS(i*N,i*N) = -2*L1(i,N)/dx^2 - 2*L2(i,N)/dzeta^2 - ...
                    L3(i,N)*Ls(i+1,N)/(4*dx*dzeta) - ...
                    L3(i,N)*Ls(i-1,N)/(4*dx*dzeta);
    LHS(i*N,(i-1)*N) = L1(i,N)/dx^2 - L2(i,N)*Ls(i,N)/dzeta^2 - ...
                    L4(i,N)*Ls(i,N)/(2*dzeta) - L5(i,N)/(2*dx);
    LHS(i*N,i*N-1) = 2*L2(i,N)/dzeta^2;
    if i<M_s-1
        LHS(i*N,(i+2)*N) = L3(i,N)*Ls(i+1,N)/(4*dx*dzeta);     
    end
    if i>2
        LHS(i*N,(i-2)*N) = L3(i,N)*Ls(i-1,N)/(4*dx*dzeta);
    end
end
end

% 2) lateral side: when i=1&M_s, u_s(i,j) = 0
%u_bdry = [0, 8.4370, 12.7761, 14.2722, 14.4966, 14.3570]; %rec
%u_bdry = [0, 6.0711, 9.1083, 10.0869, 10.2052, 10.0995]; %U convex
u_bdry = [0, 14.1157, 20.8160, 23.3082, 23.8838, 23.9179]; %U concave
u_bdry_i = interp1(0:1/5:1,u_bdry,0:1/(N-1):1);

for j = 1:N
    LHS(j,j) = 1;
    LHS((M_s-1)*N+j,(M_s-1)*N+j) = 1;
    %RHS(j,1) = u_bdry_i(j);
    %RHS((M_s-1)*N+j,1) = u_bdry_i(j);
    RHS(j,1) = 0;
    RHS((M_s-1)*N+j,1) = 0;
end

% 3) when j=1, u_s(i,j)=0
beta = 1e4;
for i = 2:M_s-1
    if abs(T_s(i,1)-(273.15-H_s(i)*8.7e-4))<1e-5
        LHS((i-1)*N+1,(i-1)*N+1+N) = L1(i,1)/dx^2+L5(i,1)/(2*dx)+L3(i,1)/(2*dx*dzeta)*H_s(i+1)*dzeta*beta/visc_s(i+1,1);
        LHS((i-1)*N+1,(i-1)*N+1) = -2*L1(i,1)/dx^2-2*L2(i,1)/dzeta^2-2*H_s(i)*dzeta*beta/visc_s(i,1)*(L2(i,1)/dzeta^2-L4(i,1)/(2*dzeta));
        LHS((i-1)*N+1,(i-1)*N+1-N) = L1(i,1)/dx^2-L5(i,1)/(2*dx)-L3(i,1)/(2*dx*dzeta)*H_s(i-1)*dzeta*beta/visc_s(i-1,1);
        LHS((i-1)*N+1,(i-1)*N+1+1) = 2*L2(i,1)/dzeta^2;
    else
        LHS((i-1)*N+1,(i-1)*N+1) = 1;
        RHS((i-1)*N+1,1) = 0;
    end
end

%for i = 2:M_s-1
%    LHS((i-1)*N+1,(i-1)*N+1) = 1;
%    RHS((i-1)*N+1,1) = 0;
%end
%------------------------ end boundary condition -------------------------
