tic

global rho g n kc Cp SPY ...
        M M_s N xi dx dzeta dzetadx dzetadx_s zeta hB hB_s hS H H_s dhSdx dhSdx_s dt ...
        de0 Sigma0

set_ice_parameters();
build_ice_geometry();

de2_s = zeros(M_s,N);

%T_s = 268.15 + zeros(M_s,N);
T_init = get_ice_T_init();
T_last = T_init;
T_s = [T_init(1,:);(T_init(1:end-1,:)+T_init(2:end,:))/2;T_init(end,:)];
A_s = 1.916e3*exp(-139000./(8.314*T_s))*(365*24*3600);

% TEST for constant LHS
%A_s = zeros(M_s,N)+1e-16;

u_s = zeros(M_s,N); 
    % velocity distribution
visc_s = 1e13/(365*24*3600)+zeros(M_s,N); 
    % initial Viscosity [Pa a]


PHTG = 0;

J = N;

iter = 0;
for yn = 1:1
%while 1
    iter = iter + 1;

    iter_u = 0;
    fprintf('yn = %d\n',yn);

while 1

    iter_u = iter_u + 1;
    
    [u_s, u] = get_ice_u(A_s, visc_s, T_s);
%    if iter_u == 1
%    [u_s, u] = get_ice_u(A_s, visc_s, T_s);
%    else
%    [u_s, u] = get_ice_u_periodic(A_s, visc_s, T_s, u_lst);
%    end

    
    [w_s, w] = get_ice_w(u_s, u);

    if 0

    if iter_u > 2
        us = u; % First, we prescribe a value us from Picard iteration
        Cs = us - u_lst;
        Sita = acos(Cs'*C/(norm(Cs)*norm(C)));        

        if isequal(Sita<=pi/8,ones(J,J))
            mu = 2.5;
        elseif isequal(Sita>pi/8,ones(J,J))&&isequal(Sita<19*pi/20,ones(J,J))
            mu = 1;
        elseif isequal(Sita>=19*pi/20,ones(J,J))
            mu = 0.5;
        end
        u = u_lst + mu*Cs;

        fprintf('error: %f\n', norm(us - u_lst)/norm(us));
        if norm(us - u_lst)/norm(us)<1e-4
            fprintf('error: %f\n', norm(us - u_lst)/norm(us));
            break;
        end
    end

    if iter_u > 1
        C = u - u_lst;
    end

    end
    
    if iter_u > 1
        relative_norm_error = 2*norm(u_s-u_s_last,'fro')/norm(u_s+u_s_last,'fro');
        fprintf('relative_norm_error: %f %f %f\n', norm(u_s-u_s_last), norm(u_s+u_s_last), relative_norm_error);
        if relative_norm_error < 1.0e-2 || iter_u == 1000
            fprintf('V solver converge! %d\n', iter_u);
            break;
        end
    end

    u_s_last = u_s;
    u_lst = u;

    [visc_s, heat_s] = get_viscosity(A_s, u_s);

    iter_T = 0;
    while 1
        iter_T = iter_T + 1;
        if iter_T > 1
            T_last = T;
        end
        heat = [heat_s(1,:);(heat_s(2:M_s-2,:)+heat_s(3:M_s-1,:))/2;heat_s(M_s,:)];
        [A, T, W] = get_ice_T(u, u_s, w, w_s, heat, T_last);
        if iter_T > 1
            rne_T = 2*norm(T-T_last)/norm(T+T_last);
            if rne_T < 1e-2 || iter_T == 2
                fprintf('T solver converges! %d\n', iter_T);
                break;
            end
        end
        T_s = [T(1,:);(T(1:end-1,:)+T(2:end,:))/2;T(end,:)];
        A_s = [A(1,:);(A(1:end-1,:)+A(2:end,:))/2;A(end,:)]*SPY;
        [visc_s, heat_s] = get_viscosity(A_s, u_s);
    end
end

% import FSM u data
% nz = 6;
% u_data = dlmread('/home/tong/Documents/myWork/2014/GlacModelComp2D/results/ESD/check_2D_T_model/convex/V10_T_model_2D/CFL_u.txt');
% %T_data = dlmread('2U/CFL_T.txt');
% 
% z = u_data(:,2);
% 
% CFL_x = u_data(z>0,1);
% CFL_z = u_data(z>0,2);
% CFL_u = u_data(z>0,3);
% CFL_v = u_data(z>0,4);
% CFL_w = u_data(z>0,5);
% 
% %CFL_T = T_data(z>0,3);
% 
% CFL_x_FS = reshape(CFL_x,nz,length(CFL_x)/nz);
% CFL_z_FS = reshape(CFL_z,nz,length(CFL_z)/nz);
% CFL_u_FS = reshape(CFL_u,nz,length(CFL_u)/nz);
% CFL_v_FS = reshape(CFL_v,nz,length(CFL_v)/nz);
% CFL_w_FS = reshape(CFL_w,nz,length(CFL_w)/nz);
% 
% z0 = 0:0.2:1;
% z1 = 0:1/(N-1):1;
% 
% u = zeros(N,M);
% w = zeros(N,M);
% for i = 1:M
%     u(:,i) = interp1(z0',CFL_u_FS(:,i), z1');
%     w(:,i) = interp1(z0',CFL_w_FS(:,i), z1');
% end
% 
% u_s = [u(:,1),(u(:,2:M)+u(:,1:M-1))/2,u(:,M)];
% w_s = [w(:,1),(w(:,2:M)+w(:,1:M-1))/2,w(:,M)];
% 
% u_s = u_s'; w_s = w_s';
% u = u'; w = w';
% 
% 
% [w_s, w] = get_ice_w(u_s, u);
% 
% 
% [visc_s, heat_s] = get_viscosity(A_s, u_s);

if 0
    if iter > 1
        relative_norm_error = 2*norm(u_s-u_s_last1)/norm(u_s+u_s_last1);
        fprintf('relative_norm_error: %f\n', relative_norm_error);
        if relative_norm_error < 1.0e-3 || iter == 2
            fprintf('V solver converge! %d\n', iter);
            break;
        end
    end

    u_s_last1 = u_s;

        heat = [heat_s(1,:);(heat_s(2:M_s-2,:)+heat_s(3:M_s-1,:))/2;heat_s(M_s,:)];
        T_last = [T_s(1,:);(T_s(2:end-2,:)+T_s(3:end-1,:))/2;T_s(end,:)];
        [A, T, W] = get_ice_T(u, u_s, w, w_s, heat_s, T_last);
        T_s = [T(1,:);(T(1:end-1,:)+T(2:end,:))/2;T(end,:)];
        A_s = [A(1,:);(A(1:end-1,:)+A(2:end,:))/2;A(end,:)]*SPY;
        [visc_s, heat_s] = get_viscosity(A_s, u_s);

end
    iter_T = 0;
    while 0
        iter_T = iter_T + 1;
        if iter_T > 1
            T_last = T;
        end
        heat = [heat_s(1,:);(heat_s(2:M_s-2,:)+heat_s(3:M_s-1,:))/2;heat_s(M_s,:)];
        [A, T, W] = get_ice_T(u, u_s, w, w_s, heat, T_last);
        if iter_T > 0
            rne_T = 2*norm(T-T_last)/norm(T+T_last);
            if rne_T < 1e-3 || iter_T == 1
                fprintf('T solver converges! %d\n', iter_T);
                break;
            end
        end
        T_s = [T(1,:);(T(1:end-1,:)+T(2:end,:))/2;T(end,:)];
        A_s = [A(1,:);(A(1:end-1,:)+A(2:end,:))/2;A(end,:)]*SPY;
        [visc_s, heat_s] = get_viscosity(A_s, u_s);
    end
end

toc

%T = T_last;

% 
 %for i = 1:M
 %    for j = 1:N
 %        A_T_pmp(i,j) = 273.15 - (N-j)*dzeta*H(i)*8.7e-4;
 %    end
 %end
 
 %A_T = A_T-A_T_pmp';
 %clear A_T_pmp;

if 1    
x_plot = ones(N,1)*xi(1:M)'/1000;
y_plot = zeta*H(1:M)'+ones(N,1)*hB(1:M)';
z_plot_T = T'-273.15; 
z_plot_uc = u';
%z_plot_P_Adv = P_Adv;
%z_plot_P_Cond = P_Cond;
%z_plot_P_Str = P_Str;
figure(1)
contourf(x_plot,y_plot,z_plot_T,min(min(T))-273.15:0.5:max(max(T))-273.15)
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
contourf(x_plot,y_plot,z_plot_uc,min(min(u)):2:max(max(u)))
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
