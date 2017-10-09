function [w_s, w] = get_ice_w(u_s, u)

global rho g n kc Cp SPY ...
        M M_s N xi dx dzeta dzetadx dzetadx_s zeta hB hB_s hS H H_s dhSdx dhSdx_s dt ...
        de0 Sigma0

U = zeros(M, N);
w_s = zeros(M, N);
w = zeros(M, N);

for i = 2:M-1
    for j = 2:N-1

        U(i,j) = (u_s(i+1,j)-u_s(i,j))/dx + dzetadx(i,j)*(u(i,j+1)-u(i,j-1))/(2*dzeta);
        U(1,j) = (u_s(1+1,j)-u_s(1,j))/(dx/2) + dzetadx(1,j)*(u(1,j+1)-u(1,j-1))/(2*dzeta);
        U(M,j) = (u_s(M+1,j)-u_s(M,j))/(dx/2) + dzetadx(M,j)*(u(M,j+1)-u(M,j-1))/(2*dzeta);
        U(i,1) = (u_s(i+1,1)-u_s(i,1))/dx + dzetadx(i,1)*(u(i,1+1)-u(i,1))/(dzeta);
        U(1,1) = (u_s(1+1,1)-u_s(1,1))/(dx/2) + dzetadx(1,1)*(u(1,1+1)-u(1,1))/(dzeta);
        U(M,1) = (u_s(M+1,1)-u_s(M,1))/(dx/2) + dzetadx(M,1)*(u(M,1+1)-u(M,1))/(dzeta);
        U(i,N) = (u_s(i+1,N)-u_s(i,N))/dx + dzetadx(i,N)*(u(i,N)-u(i,N-1))/(dzeta);
        U(1,N) = (u_s(1+1,N)-u_s(1,N))/(dx/2) + dzetadx(1,N)*(u(1,N)-u(1,N-1))/(dzeta);
        U(M,N) = (u_s(M+1,N)-u_s(M,N))/(dx/2) + dzetadx(M,N)*(u(M,N)-u(M,N-1))/(dzeta);
        
%        U(i,j) = (u(i+1,j)-u(i-1,j))/(2*dx) + dzetadx(i,j)*(u(i,j+1)-u(i,j-1))/(2*dzeta);
%        U(1,j) = (u(1+1,j)-u(1,j))/(dx) + dzetadx(1,j)*(u(1,j+1)-u(1,j-1))/(2*dzeta);
%        U(M,j) = (u(M,j)-u(M-1,j))/(dx) + dzetadx(M,j)*(u(M,j+1)-u(M,j-1))/(2*dzeta);
%        U(i,1) = (u(i+1,1)-u(i-1,1))/(2*dx) + dzetadx(i,1)*(u(i,1+1)-u(i,1))/(dzeta);
%        U(1,1) = (u(1+1,1)-u(1,1))/(dx) + dzetadx(1,1)*(u(1,1+1)-u(1,1))/(dzeta);
%        U(M,1) = (u(M,1)-u(M-1,1))/(dx) + dzetadx(M,1)*(u(M,1+1)-u(M,1))/(dzeta);
%        U(i,N) = (u(i+1,N)-u(i-1,N))/(2*dx) + dzetadx(i,N)*(u(i,N)-u(i,N-1))/(dzeta);
%        U(1,N) = (u(1+1,N)-u(1,N))/(dx) + dzetadx(1,N)*(u(1,N)-u(1,N-1))/(dzeta);
%        U(M,N) = (u(M,N)-u(M-1,N))/(dx) + dzetadx(M,N)*(u(M,N)-u(M,N-1))/(dzeta);
%

    end
end

for i = 2:M-1
    for j = 2:N

        w_s(1,1) = u(1,1)*(hB(1+1)-hB(1))/(dx) - H(1)*(0.5*U(1,1) + 0*sum(U(1,1:2)))*dzeta;
        w_s(M,1) = u(M,1)*(hB(M)-hB(M-1))/(dx) - H(M)*(0.5*U(M,1) + 0*sum(U(M,1:2)))*dzeta;
        w_s(i,1) = u(i,1)*(hB(i+1)-hB(i-1))/(2*dx) - H(i)*(0.5*U(i,1) + 1*sum(U(i,1:2)))*dzeta;
        w_s(1,j) = u(1,1)*(hB(1+1)-hB(1))/(dx) - H(1)*(0.5*U(1,1) + sum(U(1,2:j)))*dzeta;
        w_s(M,j) = u(M,1)*(hB(M)-hB(M-1))/(dx) - H(M)*(0.5*U(M,1) + sum(U(M,2:j)))*dzeta;
        w_s(i,j) = u(i,1)*(hB(i+1)-hB(i-1))/(2*dx) - H(i)*(0.5*U(i,1) + sum(U(i,2:j)))*dzeta;

    end
end

for i = 2:M-1

    w(1,:) = [u(1,1)*(hB(1+1)-hB(1))/(dx), (w_s(1,1:end-1)+w_s(1,2:end))/2];
    w(M,:) = [u(M,1)*(hB(M)-hB(M-1))/(dx), (w_s(M,1:end-1)+w_s(M,2:end))/2];
    w(i,:) = [u(i,1)*(hB(i+1)-hB(i-1))/(2*dx), (w_s(i,1:end-1)+w_s(i,2:end))/2];

end

