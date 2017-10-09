A_u = get_velocity()

R = zeros(M*N,1);
dzetadx = zeros(M,N);
L1 = zeros(M,N);
L2 = zeros(M,N);
L3 = zeros(M,N);
L4 = zeros(M,N);
L5 = zeros(M,N);
Ls = zeros(M,N);
de2 = zeros(M,N);

T_ini = 263.16 + zeros(M,N);
A0 = 1.916e3*exp(-139000./(8.314*T_ini))*(365*24*3600);

% TEST for constant A
A0 = zeros(M,N)+1e-16;

u = zeros(M,N); 
    % velocity distribution
Visc = 1e13/(365*24*3600)+zeros(M,N); 
    % initial Viscosity [Pa a]
A = zeros(M*N, M*N);

loop = 0;CALCULATE_U = 1;

%==================== Right Hand Side Matrix =========================
for i = 1:M
    for j = 1:N
      
        k = (i-1)*N + j;
        f = [0, 0.313, 0.558, 0.790, 0.884, 0.929, 0.954, 0.990, 1];
        ar = [0, 0.5, 1, 2, 3, 4, 5, 10, 1e4];
        fi = interp1(ar, f, 800/H(i));
        R(k,1) = rho*g*dhSdx(i);
            
    end
end
%=====================================================================

while CALCULATE_U == 1
    loop = loop + 1;
    
    for i = 2:M-1
    for j = 2:N-1

        k = (i-1)*N + j; % mapping method


        
        if i == 2
        dzetadx(i,j) = -1/H(i)*(1/(3*dx)*(-4*hB(i-1)+3*hB(i)+hB(i+1)) +...
                        (j-1)*dzeta/(3*dx)*(-4*H(i-1)+3*H(i)+H(i+1)));
        dzetadx(i,N) = -1/H(i)*(1/(3*dx)*(-4*hB(i-1)+3*hB(i)+hB(i+1)) +...
                        (N-1)*dzeta/(3*dx)*(-4*H(i-1)+3*H(i)+H(i+1)));
        dzetadx(i,1) = -1/H(i)*(1/(3*dx)*(-4*hB(i-1)+3*hB(i)+hB(i+1)) +...
                        (1-1)*dzeta/(3*dx)*(-4*H(i-1)+3*H(i)+H(i+1)));
        elseif i == M-1
        dzetadx(i,j) = -1/H(i)*(1/(3*dx)*(-hB(i-1)-3*hB(i)+4*hB(i+1)) +...
                        (j-1)*dzeta/(3*dx)*(-H(i-1)-3*H(i)+4*H(i+1)));
        dzetadx(i,N) = -1/H(i)*(1/(3*dx)*(-hB(i-1)-3*hB(i)+4*hB(i+1)) +...
                        (j-1)*dzeta/(3*dx)*(-H(i-1)-3*H(i)+4*H(i+1)));
        dzetadx(i,1) = -1/H(i)*(1/(3*dx)*(-hB(i-1)-3*hB(i)+4*hB(i+1)) +...
                        (j-1)*dzeta/(3*dx)*(-H(i-1)-3*H(i)+4*H(i+1)));
        else
        dzetadx(i,j) = -((hB(i+1)-hB(i-1))/(2*dx)+...
                        (j-1)*dzeta*(H(i+1)-H(i-1))/(2*dx))/H(i);
        dzetadx(i,N) = -((hB(i+1)-hB(i-1))/(2*dx)+...
                        (N-1)*dzeta*(H(i+1)-H(i-1))/(2*dx))/H(i);
        dzetadx(i,1) = -((hB(i+1)-hB(i-1))/(2*dx)+...
                        (1-1)*dzeta*(H(i+1)-H(i-1))/(2*dx))/H(i);
        end
        dzetadx(M,N) = -((hB(M)-hB(M-1))/(dx)+...
                        (N-1)*dzeta*(H(M)-H(M-1))/(dx))/H(M);

        L1(i,j) = 4*Visc(i,j); 
        L1(i,N) = 4*Visc(i,N);
        L1(i,1) = 4*Visc(i,1);
        L1(1,j) = 4*Visc(1,j);
        

        L2(i,j) = 4*Visc(i,j)*((dzetadx(i,j))^2 + 1/4*(1/H(i))^2);
        L2(i,N) = 4*Visc(i,N)*((dzetadx(i,N))^2 + 1/4*(1/H(i))^2);
        L2(i,1) = 4*Visc(i,1)*((dzetadx(i,1))^2 + 1/4*(1/H(i))^2);
        L2(1,j) = 4*Visc(1,j)*((dzetadx(1,j))^2 + 1/4*(1/H(i))^2);

        L3(i,j) = 8*Visc(i,j)*dzetadx(i,j);
        L3(i,N) = 8*Visc(i,N)*dzetadx(i,N);
        L3(i,1) = 8*Visc(i,1)*dzetadx(i,1);
        L3(1,j) = 8*Visc(1,j)*dzetadx(1,j);

        if i == 2
        L4(i,j) = 4*Visc(i,j)*(1/(3*dx)*(-4*dzetadx(i-1,j)+3*dzetadx(i,j)+dzetadx(i+1,j)))+...
                4*Visc(i,j)*dzetadx(i,j)*(dzetadx(i,j+1)-dzetadx(i,j-1))/(2*dzeta) +...
                4*dzetadx(i,j)*(1/(3*dx)*(-4*Visc(i-1,j)+3*Visc(i,j)+Visc(i+1,j))) +...
                (Visc(i,j+1)-Visc(i,j-1))/(2*dzeta)*(4*dzetadx(i,j)^2+(1/H(i))^2);
        L4(i,N) = 4*Visc(i,N)*(1/(3*dx)*(-4*dzetadx(i-1,N)+3*dzetadx(i,N)+dzetadx(i+1,N)))+...
                4*Visc(i,N)*dzetadx(i,N)*(dzetadx(i,N)-dzetadx(i,N-1))/(dzeta) +...
                4*dzetadx(i,N)*(1/(3*dx)*(-4*Visc(i-1,N)+3*Visc(i,N)+Visc(i+1,N))) +...
                (Visc(i,N)-Visc(i,N-1))/(dzeta)*(4*dzetadx(i,N)^2+(1/H(i))^2);
        L4(i,1) = 4*Visc(i,1)*(1/(3*dx)*(-4*dzetadx(i-1,1)+3*dzetadx(i,1)+dzetadx(i+1,1)))+...
                4*Visc(i,1)*dzetadx(i,1)*(dzetadx(i,2)-dzetadx(i,1))/(dzeta) +...
                4*dzetadx(i,1)*(1/(3*dx)*(-4*Visc(i-1,1)+3*Visc(i,1)+Visc(i+1,1))) +...
                (Visc(i,2)-Visc(i,1))/(dzeta)*(4*dzetadx(i,1)^2+(1/H(i))^2);
        elseif i == M-1
        L4(i,j) = 4*Visc(i,j)*(1/(3*dx)*(-dzetadx(i-1,j)-3*dzetadx(i,j)+4*dzetadx(i+1,j)))+...
                4*Visc(i,j)*dzetadx(i,j)*(dzetadx(i,j+1)-dzetadx(i,j-1))/(2*dzeta) +...
                4*dzetadx(i,j)*(1/(3*dx)*(-Visc(i-1,j)-3*Visc(i,j)+4*Visc(i+1,j))) +...
                (Visc(i,j+1)-Visc(i,j-1))/(2*dzeta)*(4*dzetadx(i,j)^2+(1/H(i))^2);
        L4(i,N) = 4*Visc(i,N)*(1/(3*dx)*(-dzetadx(i-1,N)-3*dzetadx(i,N)+4*dzetadx(i+1,N)))+...
                4*Visc(i,N)*dzetadx(i,N)*(dzetadx(i,N)-dzetadx(i,N-1))/(dzeta) +...
                4*dzetadx(i,N)*(1/(3*dx)*(-Visc(i-1,N)-3*Visc(i,N)+4*Visc(i+1,N))) +...
                (Visc(i,N)-Visc(i,N-1))/(dzeta)*(4*dzetadx(i,N)^2+(1/H(i))^2);
        L4(i,1) = 4*Visc(i,1)*(1/(3*dx)*(-dzetadx(i-1,1)-3*dzetadx(i,1)+4*dzetadx(i+1,1)))+...
                4*Visc(i,1)*dzetadx(i,1)*(dzetadx(i,1+1)-dzetadx(i,1))/(dzeta) +...
                4*dzetadx(i,1)*(1/(3*dx)*(-Visc(i-1,1)-3*Visc(i,1)+4*Visc(i+1,1))) +...
                (Visc(i,1+1)-Visc(i,1))/(dzeta)*(4*dzetadx(i,1)^2+(1/H(i))^2);
        else
        L4(i,j) = 4*Visc(i,j)*(dzetadx(i+1,j)-dzetadx(i-1,j))/...
                    (2*dx)+ 4*Visc(i,j)*dzetadx(i,j)*...
                    (dzetadx(i,j+1)-dzetadx(i,j-1))/(2*dzeta) +...
                  4*(Visc(i+1,j)-Visc(i-1,j))/(2*dx)*dzetadx(i,j) +...
                   (Visc(i,j+1)-Visc(i,j-1))/(2*dzeta)*(4*dzetadx(i,j)^2+(1/H(i))^2);
        L4(i,N) = 4*Visc(i,N)*(dzetadx(i+1,N)-dzetadx(i-1,N))/...
                    (2*dx)+ 4*Visc(i,N)*dzetadx(i,N)*...
                    (dzetadx(i,N)-dzetadx(i,N-1))/(dzeta) + ...
                  4*(Visc(i+1,N)-Visc(i-1,N))/(2*dx)*dzetadx(i,N) +...
                   (Visc(i,N)-Visc(i,N-1))/(2*dzeta)*(4*dzetadx(i,N)^2+(1/H(i))^2);
        L4(i,1) = 4*Visc(i,1)*(dzetadx(i+1,1)-dzetadx(i-1,1))/(2*dx)+...
                     4*Visc(i,1)*dzetadx(i,1)*(dzetadx(i,2)-dzetadx(i,1))/(dzeta) + ...
                  4*(Visc(i+1,1)-Visc(i-1,1))/(2*dx)*dzetadx(i,1) +...
                   (Visc(i,2)-Visc(i,1))/(dzeta)*(4*dzetadx(i,1)^2+(1/H(i))^2);        
        L4(1,j) = 4*Visc(1,j)*(dzetadx(2,j)-dzetadx(1,j))/...
                    (dx)+ 4*Visc(1,j)*dzetadx(1,j)*...
                    (dzetadx(1,j+1)-dzetadx(1,j-1))/(2*dzeta) +...
                  4*(Visc(2,j)-Visc(1,j))/(dx)*dzetadx(1,j) +...
                   (Visc(1,j+1)-Visc(1,j-1))/(2*dzeta)*(4*dzetadx(1,j)^2+(1/H(1))^2);
        end

        if i == 2
        L5(i,j) = 4*(1/(3*dx)*(-4*Visc(i-1,j)+3*Visc(i,j)+Visc(i+1,j))) +...
                4*(Visc(i,j+1)-Visc(i,j-1))/(2*dzeta)*(dzetadx(i,j));
        L5(i,N) = 4*(1/(3*dx)*(-4*Visc(i-1,N)+3*Visc(i,N)+Visc(i+1,N))) +...
                4*(Visc(i,N)-Visc(i,N-1))/(dzeta)*(dzetadx(i,N));
        L5(i,1) = 4*(1/(3*dx)*(-4*Visc(i-1,1)+3*Visc(i,j)+Visc(i+1,1))) +...
                4*(Visc(i,2)-Visc(i,1))/(dzeta)*(dzetadx(i,1));
        elseif i == M-1
        L5(i,j) = 4*(1/(3*dx)*(-Visc(i-1,j)-3*Visc(i,j)+4*Visc(i+1,j))) +...
                4*(Visc(i,j+1)-Visc(i,j-1))/(2*dzeta)*(dzetadx(i,j));
        L5(i,N) = 4*(1/(3*dx)*(-Visc(i-1,N)-3*Visc(i,N)+4*Visc(i+1,N))) +...
                4*(Visc(i,N)-Visc(i,N-1))/(dzeta)*(dzetadx(i,N));
        L5(i,1) = 4*(1/(3*dx)*(-Visc(i-1,1)-3*Visc(i,1)+4*Visc(i+1,1))) +...
                4*(Visc(i,1+1)-Visc(i,1))/(dzeta)*(dzetadx(i,1));
        else
        L5(i,j) = 4*(Visc(i+1,j)-Visc(i-1,j))/(2*dx) + ...
            4*(Visc(i,j+1)-Visc(i,j-1))/(2*dzeta)*(dzetadx(i,j));
        L5(i,N) = 4*(Visc(i+1,N)-Visc(i-1,N))/(2*dx) + ...
            4*(Visc(i,N)-Visc(i,N-1))/(dzeta)*(dzetadx(i,N));
        L5(i,1) = 4*(Visc(i+1,1)-Visc(i-1,1))/(2*dx) + ...
            4*(Visc(i,2)-Visc(i,1))/(dzeta)*(dzetadx(i,1));
        end

        Ls(i,N) = 4*dhSdx(i)*dzeta/dx/(1/H(i)-4*dhSdx(i)*dzetadx(i,N));
        Ls(1,N) = 4*dhSdx(1)*dzeta/dx/(1/H(1)-4*dhSdx(1)*dzetadx(1,N));
        Ls(M,N) = 4*dhSdx(M)*dzeta/dx/(1/H(M)-4*dhSdx(M)*dzetadx(M,N));

        if i == 2
            A(k,k-N) = 8*L1(i,j)/(3*dx^2) - 4*L5(i,j)/(3*dx);
            A(k,k) = -4*L1(i,j)/dx^2 - 2*L2(i,j)/dzeta^2 + L5(i,j)/dx;
            A(k,k+N) = 4*L1(i,j)/(3*dx^2) + L5(i,j)/(3*dx);
            A(k,k+1) = L2(i,j)/dzeta^2 + L3(i,j)/(2*dx*dzeta) + L4(i,j)/(2*dzeta);
            A(k,k-1) = L2(i,j)/dzeta^2 - L3(i,j)/(2*dx*dzeta) - L4(i,j)/(2*dzeta);
            A(k,k-N+1) = -2*L3(i,j)/(3*dx*dzeta);
            A(k,k-N-1) = 2*L3(i,j)/(3*dx*dzeta);
            A(k,k+N+1) = L3(i,j)/(6*dx*dzeta);
            A(k,k+N-1) = -L3(i,j)/(6*dx*dzeta);
        elseif i == M-1
            A(k,k-N) = 4*L1(i,j)/(3*dx^2) - L5(i,j)/(3*dx);
            A(k,k) = -4*L1(i,j)/dx^2 - 2*L2(i,j)/dzeta^2 - L5(i,j)/dx;
            A(k,k+N) = 8*L1(i,j)/(3*dx^2) + 4*L5(i,j)/(3*dx);
            A(k,k+1) = L2(i,j)/dzeta^2 - L3(i,j)/(2*dx*dzeta) + L4(i,j)/(2*dzeta);
            A(k,k-1) = L2(i,j)/dzeta^2 + L3(i,j)/(2*dx*dzeta) - L4(i,j)/(2*dzeta);
            A(k,k-N+1) = -L3(i,j)/(6*dx*dzeta);
            A(k,k-N-1) = L3(i,j)/(6*dx*dzeta);
            A(k,k+N+1) = 2*L3(i,j)/(3*dx*dzeta);
            A(k,k+N-1) = -2*L3(i,j)/(3*dx*dzeta);
        else
            A(k,k-N-1) = L3(i,j)/(4*dx*dzeta);
            A(k,k-N) = L1(i,j)/dx^2 - L5(i,j)/(2*dx);
            A(k,k-N+1) = -L3(i,j)/(4*dx*dzeta);
            A(k,k-1) = L2(i,j)/dzeta^2 - L4(i,j)/(2*dzeta);
            A(k,k) = -2*L1(i,j)/dx^2 - 2*L2(i,j)/dzeta^2;
            A(k,k+1) = L2(i,j)/dzeta^2 + L4(i,j)/(2*dzeta);
            A(k,k+N-1) = -L3(i,j)/(4*dx*dzeta);
            A(k,k+N) = L1(i,j)/dx^2 + L5(i,j)/(2*dx);
            A(k,k+N+1) = L3(i,j)/(4*dx*dzeta);
        end

    end
    end

    %=============== boundary condition =================

    % surface: when i=2:M-1, j=N
    for i = 2:M-1
        if i == 2
            A(i*N,(i-1)*N) = 8*L1(i,N)/(3*dx^2) - 4*L5(i,N)/(3*dx) -...
                            8*Ls(i,N)/3*(L2(i,N)/dzeta^2 + L3(i,N)/(2*dx*dzeta) + L4(i,N)/(2*dzeta)) +...
                            4*L3(i,N)/(3*dx*dzeta)*Ls(i-1,N);
            A(i*N,i*N) = -4*L1(i,N)/dx^2 - 2*L2(i,N)/dzeta^2 + L5(i,N)/dx +...
                        2*Ls(i,N)*(L2(i,N)/dzeta^2 + L3(i,N)/(2*dx*dzeta) + L4(i,N)/(2*dzeta)) - ...
                        4*L3(i,N)/(9*dx*dzeta)*Ls(i-1,N) - 4*L3(i,N)/(9*dx*dzeta)*Ls(i+1,N);
            A(i*N,(i+1)*N) = 4*L1(i,N)/(3*dx^2) + L5(i,N)/(3*dx) +...
                            2*Ls(i,N)/3*(L2(i,N)/dzeta^2 + L3(i,N)/(2*dx*dzeta) + L4(i,N)/(2*dzeta)) +...
                            L3(i,N)/(3*dx*dzeta)*Ls(i+1,N);
            A(i*N,i*N-1) = L2(i,N)/dzeta^2 - L3(i,N)/(2*dx*dzeta) - L4(i,N)/(2*dzeta) +...
                            (L2(i,N)/dzeta^2 + L3(i,N)/(2*dx*dzeta) + L4(i,N)/(2*dzeta));
            A(i*N,(i+2)*N) = L3(i,N)/(9*dx*dzeta)*Ls(i+1,N);
        elseif i == M-1
            A(i*N,(i-1)*N) = 4*L1(i,N)/(3*dx^2) - L5(i,N)/(3*dx) -...
                            2*Ls(i,N)/3*(L2(i,N)/dzeta^2 - L3(i,N)/(2*dx*dzeta) + L4(i,N)/(2*dzeta)) +...
                            L3(i,N)/(3*dx*dzeta)*Ls(i-1,N);
            A(i*N,i*N) = -4*L1(i,N)/dx^2 - 2*L2(i,N)/dzeta^2 - L5(i,N)/dx -...
                        2*Ls(i,N)*(L2(i,N)/dzeta^2 - L3(i,N)/(2*dx*dzeta) + L4(i,N)/(2*dzeta)) - ...
                        4*L3(i,N)/(9*dx*dzeta)*Ls(i-1,N) - 4*L3(i,N)/(9*dx*dzeta)*Ls(i+1,N);
            A(i*N,(i+1)*N) = 8*L1(i,N)/(3*dx^2) + 4*L5(i,N)/(3*dx) -...
                            8*Ls(i,N)/3*(L2(i,N)/dzeta^2 - L3(i,N)/(2*dx*dzeta) + L4(i,N)/(2*dzeta)) -...
                            4*L3(i,N)/(3*dx*dzeta)*Ls(i+1,N);
            A(i*N,i*N-1) = L2(i,N)/dzeta^2 + L3(i,N)/(2*dx*dzeta) - L4(i,N)/(2*dzeta) +...
                            L2(i,N)/dzeta^2 - L3(i,N)/(2*dx*dzeta) + L4(i,N)/(2*dzeta);
            A(i*N,(i-2)*N) = L3(i,N)/(9*dx*dzeta)*Ls(i-1,N);
        else
            A(i*N,(i+1)*N) = L1(i,N)/dx^2 + L2(i,N)*Ls(i,N)/dzeta^2 +...
                             L4(i,N)*Ls(i,N)/(2*dzeta) + L5(i,N)/(2*dx);
            A(i*N,i*N) = -2*L1(i,N)/dx^2 - 2*L2(i,N)/dzeta^2 - ...
                            L3(i,N)*Ls(i+1,N)/(4*dx*dzeta) - ...
                            L3(i,N)*Ls(i-1,N)/(4*dx*dzeta);
            A(i*N,(i-1)*N) = L1(i,N)/dx^2 - L2(i,N)*Ls(i,N)/dzeta^2 - ...
                            L4(i,N)*Ls(i,N)/(2*dzeta) - L5(i,N)/(2*dx);
            A(i*N,i*N-1) = 2*L2(i,N)/dzeta^2;
            if i<M-1
                A(i*N,(i+2)*N) = L3(i,N)*Ls(i+1,N)/(4*dx*dzeta);     
            end
            if i>2
                A(i*N,(i-2)*N) = L3(i,N)*Ls(i-1,N)/(4*dx*dzeta);
            end
        end
    end

    % lateral side: when i=1&M, u(i,j) = 0
    for j = 1:N
        A(j,j) = 1;
        A((M-1)*N+j,(M-1)*N+j) = 1;
        R(j,1) = 0;
        R((M-1)*N+j,1) = 0;
    end

    % when j=1, u(i,j)=0
    for i = 1:M
        A((i-1)*N+1,(i-1)*N+1) = 1;
        R((i-1)*N+1,1) = 0;
    end
    %==========================================================

    v = A\R;

    for k=1:M*N
        if fix(k/N) == k/N
            i = k/N;
        else
            i = fix(k/N) + 1;
        end
        j = k-(i-1)*N;
        u(i,j) = v(k);
    end

    u(1,:) = 0; u(M,:) = 0; u(:,1) = 0;
    loop
    if loop == 1; CALCULATE_U = 0; end
        
    if loop > 1
        relative_norm_error = 2*norm(u-u_lst)/norm(u+u_lst);
        if relative_norm_error < 1e-3
            fprintf('V solver converge! %d\n', loop);
            CALCULATE_U = 0;
        end
    end

    u_lst = u;

    u_o = [u(1,:);(u(2:M-2,:)+u(3:M-1))/2;u(M,:)];

    for i = 1:M
        for j = 1:N

            im1 = i-1; if i==1; im1 = 1;end
            ip1 = i+1; if i==M; ip1 = M;end
            jm1 = j-1; if j==1; jm1 = 1;end
            jp1 = j+1; if j==N; jp1 = N;end
            dzeta2 = 2*dzeta; if j==1||j==N; dzeta2 = dzeta; end
            dx2 = 2*dx; if i==1||i==M; dx2 = dx; end

            if i == 2
            de2(i,j)=1/4/H(i)^2*((u(i,jp1)-u(i,jm1))/(dzeta2))^2+... 
            	(1/(3*dx)*(-4*u(im1,j)+3*u(i,j)+u(ip1,j)) + ...
            	dzetadx(i,j)*(u(i,jp1)-u(i,jm1))/(dzeta2))^2;
            elseif j == M-1
            de2(i,j)=1/4/H(i)^2*((u(i,jp1)-u(i,jm1))/(dzeta2))^2+... 
            	(1/(3*dx)*(-u(im1,j)-3*u(i,j)+4*u(ip1,j)) + ...
            	dzetadx(i,j)*(u(i,jp1)-u(i,jm1))/(dzeta2))^2;
            else
            de2(i,j)=1/4/H(i)^2*((u(i,jp1)-u(i,jm1))/(dzeta2))^2+... 
            	((u(ip1,j)-u(im1,j))/dx2 + ...
            	dzetadx(i,j)*(u(i,jp1)-u(i,jm1))/(dzeta2))^2;

            Visc(i,j)=1/2*A0(i,j)^(-1/n)*(de2(i,j)+de0)^((1-n)/(2*n));

        end
    end
    
    %clear A0;
    
    A_uc = u';A_Viscc = Visc';A_de2c = de2'; T = T_ini';

    iter_T = 1;
    %[A0, A_T, A_w, h_CTS]=Temperature_2D(A_uc, A_Viscc, A_de2c, T);
    while 0
    
        T_last = T;
        [A0, A_T, A_w, h_CTS]=Temperature_2D(A_uc, A_Viscc, A_de2c, T);
        T = A_T;

        if iter_T > 1
            relative_norm_error_T=2*norm(A_T-T_last)/norm(A_T+T_last);
            if relative_norm_error_T < 1e-8
                fprintf('T solver converge! %d\n', iter_T);
                break;
            end
        end

        iter_T = iter_T+1;
    end

    
    %A0 = A0'*(365*24*3600);
    %T_ini = A_T'; A_w = A_w * (365*24*3600);
    
    for i = 1:M
        for j = 1:N

            im1 = i-1; if i==1; im1 = 1;end
            ip1 = i+1; if i==M; ip1 = M;end
            jm1 = j-1; if j==1; jm1 = 1;end
            jp1 = j+1; if j==N; jp1 = N;end
            dzeta2 = 2*dzeta; if j==1||j==N; dzeta2 = dzeta; end
            dx2 = 2*dx; if i==1||i==M; dx2 = dx; end

            de2(i,j)=1/4/H(i)^2*((u(i,jp1)-u(i,jm1))/(dzeta2))^2 +...
            		 ((u(ip1,j)-u(im1,j))/(dx2) + ...
            		 dzetadx(i,j)*(u(i,jp1)-u(i,jm1))/(dzeta2))^2;

            Visc(i,j)=1/2*A0(i,j)^(-1/n)*(de2(i,j)+de0)^((1-n)/(2*n));

        end
    end
    
     
end

    %[A0, A_T, A_w, h_CTS]=Temperature_2D(A_uc, A_Viscc, A_de2c, T);
while 0

    for i = 1:M
        for j = 1:N

            im1 = i-1; if i==1; im1 = 1;end
            ip1 = i+1; if i==M; ip1 = M;end
            jm1 = j-1; if j==1; jm1 = 1;end
            jp1 = j+1; if j==N; jp1 = N;end
            dzeta2 = 2*dzeta; if j==1||j==N; dzeta2 = dzeta; end
            dx2 = 2*dx; if i==1||i==M; dx2 = dx; end

            de2(i,j)=1/4/H(i)^2*((u(i,jp1)-u(i,jm1))/(dzeta2))^2 +...
            		 ((u(ip1,j)-u(im1,j))/(dx2) + ...
            		 dzetadx(i,j)*(u(i,jp1)-u(i,jm1))/(dzeta2))^2;

            Visc(i,j)=1/2*A0(i,j)^(-1/n)*(de2(i,j)+de0)^((1-n)/(2*n));

        end
    end
    A_Viscc = Visc';

    T_last = T;
    [A0, A_T, A_w, h_CTS]=Temperature_2D(A_uc, A_Viscc, A_de2c, T);
    A0 = A0'*(365*24*3600);
    T = A_T;

    if iter_T > 1
        relative_norm_error_T=2*norm(A_T-T_last)/norm(A_T+T_last);
        if relative_norm_error_T < 1e-16
            fprintf('T solver converge! %d\n', iter_T);
            break;
        end
    end

    iter_T = iter_T+1;
end

% after V converges, we update hSc

%hSc_last = hSc;
%hSc = Thickness_2D(hSc_last,hBc,A_uc);

%if 2*norm(hSc-hSc_last)/norm(hSc+hSc_last) < 1e-4
%    fprintf('geometry converge!\n');
%    break;
%end

end     

toc

% 
 for i = 1:M
     for j = 1:N
         A_T_pmp(i,j) = 273.15 - (N-j)*dzeta*H(i)*8.7e-4;
     end
 end
 
 A_T = A_T-A_T_pmp';
 clear A_T_pmp;

if 0    
x_plot = ones(N,1)*xic(1:M)'/1000;
y_plot = zeta*H(1:M)'+ones(N,1)*hBc(1:M)';
z_plot_T = T(:,1:M)-273.15; 
z_plot_uc = A_uc(:,1:M);
%z_plot_P_Adv = P_Adv;
%z_plot_P_Cond = P_Cond;
%z_plot_P_Str = P_Str;
figure(1)
contourf(x_plot,y_plot,z_plot_T,min(min(T))-273.15:5:max(max(T))-273.15)
%hold on
%plot(y/1000, h_CTS);
xlabel('Glacier length [km]','Fontsize',30)
ylabel('Altitude [m]','Fontsize',30)
colorbar('EastOutside')                
%hold on
%fill([y'/1000 y(end)/1000 y(1)/1000],[hBc' 5300 5300],'k');
set(gca,'FontSize',30)
%plot(x(1:I)./1000,h_CTS(1:I),'o','LineWidth',2)
figure(2)
contourf(x_plot,y_plot,z_plot_uc,min(min(A_uc)):100:max(max(A_uc)))
xlabel('Horiztontal distance (km)','Fontsize',30)
ylabel('Elevation (m)','Fontsize',30)
%axis([0 x(end)/1000 5300 6600])
%hold on
%fill([y'/1000 y(end)/1000 y(1)/1000],[hBc' 5300 5300],'k');
colorbar('EastOutside')
set(gca,'FontSize',30) 

%figure(3)
%plot(hSc,A_uc(end,:),'--k')
%hold on
%plot(A_T(:,30)-273.15,H(30):-dzeta*H(30):0)
%hold on
%plot(A_T(:,41)-273.15,H(41):-dzeta*H(41):0)
%hold on
%plot(A_T(:,43)-273.15,H(43):-dzeta*H(43):0)
%hold on
%plot(A_T(:,45)-273.15,H(45):-dzeta*H(45):0)
%hold on
%plot(A_T(:,51)-273.15,H(51):-dzeta*H(51):0)
%hold on
%plot(-7.24,17.5,'o',-1.86,18,'o','Markersize',16)

%figure(5)
%contourf(x_plot,y_plot,Pe,min(min(Pe)):5:max(max(Pe)))
%xlabel('Glacier length [km]')
%ylabel('Altitude [m]')
%colorbar('EastOutside')                
%hold on
%fill([y'/1000 y(end)/1000 y(1)/1000],[hBc' 5300 5300],'k');
end
