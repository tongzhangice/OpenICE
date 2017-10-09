function [visc_s, heat_s] = get_viscosity(A_s, u_s)

global rho g n kc Cp SPY ...
        M M_s N xi dx dzeta dzetadx dzetadx_s zeta hB hB_s hS H H_s dhSdx dhSdx_s dt ...
        de0 Sigma0

for i = 1:M_s
for j = 1:N

    im1 = i-1; if i==1; im1 = 1;end
    ip1 = i+1; if i==M_s; ip1 = M_s;end
    jm1 = j-1; if j==1; jm1 = 1;end
    jp1 = j+1; if j==N; jp1 = N;end
    dzeta2 = 2*dzeta; if j==1||j==N; dzeta2 = dzeta; end
    dx2 = 2*dx; if i==1||i==M_s; dx2 = dx; end

    if i == 1
    de2_s(1,j)=1/4/H_s(1)^2*((u_s(1,jp1)-u_s(1,jm1))/(dzeta2))^2+... 
        ((u_s(2,j)-u_s(1,j))/(dx/2) + ...
        dzetadx_s(1,j)*(u_s(1,jp1)-u_s(1,jm1))/(dzeta2))^2;
    elseif i == 2
     de2_s(i,j)=1/4/H_s(i)^2*((u_s(i,jp1)-u_s(i,jm1))/(dzeta2))^2+... 
         (1/(3*dx)*(-4*u_s(im1,j)+3*u_s(i,j)+u_s(ip1,j)) + ...
         dzetadx_s(i,j)*(u_s(i,jp1)-u_s(i,jm1))/(dzeta2))^2;
    elseif i == M_s-1
    de2_s(i,j)=1/4/H_s(i)^2*((u_s(i,jp1)-u_s(i,jm1))/(dzeta2))^2+... 
        (1/(3*dx)*(-u_s(im1,j)-3*u_s(i,j)+4*u_s(ip1,j)) + ...
        dzetadx_s(i,j)*(u_s(i,jp1)-u_s(i,jm1))/(dzeta2))^2;
    elseif i == M_s
    de2_s(M_s,j)=1/4/H_s(M_s)^2*((u_s(M_s,jp1)-u_s(M_s,jm1))/(dzeta2))^2+... 
        ((u_s(M_s,j)-u_s(M_s-1,j))/(dx/2)+ ...
        dzetadx_s(M_s,j)*(u_s(M_s,jp1)-u_s(M_s,jm1))/(dzeta2))^2;
    else
    de2_s(i,j)=1/4/H_s(i)^2*((u_s(i,jp1)-u_s(i,jm1))/(dzeta2))^2+... 
        ((u_s(ip1,j)-u_s(im1,j))/dx2 + ...
        dzetadx_s(i,j)*(u_s(i,jp1)-u_s(i,jm1))/(dzeta2))^2;
    end

    if i == 1
    de2_s1(1,j)=1/4/H_s(1)^2*((u_s(1,jp1)-u_s(1,jm1))/(dzeta2))^2+... 
        ((u_s(2,j)-u_s(1,j))/(dx/2)*0 + ...
        dzetadx_s(1,j)*(u_s(1,jp1)-u_s(1,jm1))/(dzeta2)*0)^2;
     elseif i == 2
     de2_s1(i,j)=1/4/H_s(i)^2*((u_s(i,jp1)-u_s(i,jm1))/(dzeta2))^2+... 
         (1/(3*dx)*(-4*u_s(im1,j)+3*u_s(i,j)+u_s(ip1,j))*0 + ...
         dzetadx_s(i,j)*(u_s(i,jp1)-u_s(i,jm1))/(dzeta2)*0)^2;
   elseif i == M_s-1
    de2_s1(i,j)=1/4/H_s(i)^2*((u_s(i,jp1)-u_s(i,jm1))/(dzeta2))^2+... 
        (1/(3*dx)*(-u_s(im1,j)-3*u_s(i,j)+4*u_s(ip1,j))*0 + ...
        dzetadx_s(i,j)*(u_s(i,jp1)-u_s(i,jm1))/(dzeta2)*0)^2;
    elseif i == M_s
    de2_s1(M_s,j)=1/4/H_s(M_s)^2*((u_s(M_s,jp1)-u_s(M_s,jm1))/(dzeta2))^2+... 
        ((u_s(M_s,j)-u_s(M_s-1,j))/(dx/2)*0+ ...
        dzetadx_s(M_s,j)*(u_s(M_s,jp1)-u_s(M_s,jm1))/(dzeta2)*0)^2;
    else
    de2_s1(i,j)=1/4/H_s(i)^2*((u_s(i,jp1)-u_s(i,jm1))/(dzeta2))^2+... 
        ((u_s(ip1,j)-u_s(im1,j))/dx2*0 + ...
        dzetadx_s(i,j)*(u_s(i,jp1)-u_s(i,jm1))/(dzeta2)*0)^2;
    end
    visc_s(i,j)=1/2*A_s(i,j)^(-1/n)*(de2_s(i,j)+de0)^((1-n)/(2*n));

    heat_s(i,j) = 4*visc_s(i,j)*de2_s(i,j); 

end
end
