function build_ice_geometry()

global rho g n kc Cp SPY ...
        M M_s N xi dx dzeta dzetadx dzetadx_s zeta hB hB_s hS H H_s dhSdx dhSdx_s dt ...
        de0 Sigma0

M = 61;
    % divide x direction into M-1 intervals                                                                 
N = 21;
    % divide z direction into N-1 intervals
xi = [0:100:12100]';
%hS = 3101 - 0.2*xi;
%hB = zeros(M, 1);

%hB = hS - 100 - 15*abs(sin(2*pi/2000*xi));
% for i = 1:M
%     if xi(i) >= 1600 && xi(i) <= 2400
%         %hB(i) = hS(i) - 100 + 10*sin(2*pi/1600*(xi(i)-1200));
%         hB(i) = hS(i) - 100 + 40*exp(-(xi(i)-2000)^2/1000^2);
% %        hB(i) = hS(i) - 100;
%     else
%         hB(i) = hS(i) - 100;
%     end
% end


%hB = hS - 100;
%H = xi*0+100;

geo = dlmread('arolla100.dat');
%xi = geo(:,1);
%xi1 = [0:20:5000]';
hB = geo(:,2);
%hB = interp1(xi,hB,xi1);
hS = geo(:,3)+1.0;
%hS = interp1(xi,hS,xi1);


load hSc.mat;
load hBc.mat;

hS = hSc'; hB = hBc';

%geo = dlmread('EISMINT_F_2D_STEADY_GEO.txt');
%xi = geo(:,1)*1e3;
%H = geo(:,2)*1e3;
%xic = geo(:,1); xic1 = [0:50:5000]';
%hBc = geo(:,2); hBc = interp1(xic, hBc, xic1);
%hSc = geo(:,3)+0.1; hSc = interp1(xic, hSc, xic1);
%H = hS - hB;

%xi = [0:50e3:750e3]';
%xi = [0:1e2:6e3]';

dx = max(xi)/(M-1);
    % x spacing [m]
dzeta = 1/(N-1); zeta = 0:dzeta:1; zeta = zeta';

%y = 0:dx:max(xic); y=y';
    % Distance bedrock position [m];

%hB = 3000-0.2*xi;
%hB = xi*0; %hB = hBc;
    % Bedrock elevation [m]
%hB = interp1(xi, hB, xi1);
%H = interp1(xi, H, xi1);
%H = 2*0.02*(xi.*(8000-xi)).^0.5 + 1;
%H = xi.*0.0 + 101; 
%hS = hB + H;
%hS = 3571.1*(1-(xi/750e3).^(1+1/n)).^(1/(2+2/n))+1;
%xh = [0,100,200,300,400,500,600,700,750]*1e3';
%hh = [3522,3411.8,3264.7,3058.8,2794.1,2470.6,2000,1176.5,1]';
%EISMINT_F = dlmread('EISMINT_F_2D_STEADY_GEO.txt');
%xh = EISMINT_F(:,1)*1e3;xh=xh(1:24);
%hh = EISMINT_F(:,2)*1e3;hh=hh(1:24);
%hS = interp1(xh,hh,xi);

H = hS - hB;

dhSdx_s = zeros(M+1,1);
dhSdx_s(2:M) = (hS(2:M) - hS(1:M-1))/dx; 
dhSdx_s(1) = 0; dhSdx_s(M+1) = dhSdx_s(M);

H_s = [H(1);(H(1:M-1)+H(2:M))/2;H(M)];
hB_s = [hB(1);(hB(1:M-1)+hB(2:M))/2;hB(M)];
%hS = hB_s + H_s;

M_s = M+1;

dzetadx_s = zeros(M_s,N);
dzetadx = zeros(M,N);
PHTG = 1;

%for i = 1:M_s
%for j = 2:N
%
%    if i == 1
%    dzetadx_s(i,j) = -((hB_s(i+1)-hB_s(1))/(dx/2)+...
%                    (j-1)*dzeta*(H_s(i+1)-H_s(1))/(dx/2))/H_s(i);
%    elseif i == 2 && PHTG
%    dzetadx_s(i,j) = -1/H_s(i)*(1/(3*dx)*(-4*hB_s(i-1)+3*hB_s(i)+hB_s(i+1)) +...
%                    (j-1)*dzeta/(3*dx)*(-4*H_s(i-1)+3*H_s(i)+H_s(i+1)));
%    dzetadx_s(i,N) = -1/H_s(i)*(1/(3*dx)*(-4*hB_s(i-1)+3*hB_s(i)+hB_s(i+1)) +...
%                    (N-1)*dzeta/(3*dx)*(-4*H_s(i-1)+3*H_s(i)+H_s(i+1)));
%    dzetadx_s(i,1) = -1/H_s(i)*(1/(3*dx)*(-4*hB_s(i-1)+3*hB_s(i)+hB_s(i+1)) +...
%                    (1-1)*dzeta/(3*dx)*(-4*H_s(i-1)+3*H_s(i)+H_s(i+1)));
%    elseif i == M_s-1 && PHTG && 1
%    dzetadx_s(i,j) = -1/H_s(i)*(1/(3*dx)*(-hB_s(i-1)-3*hB_s(i)+4*hB_s(i+1)) +...
%                    (j-1)*dzeta/(3*dx)*(-H_s(i-1)-3*H_s(i)+4*H_s(i+1)));
%    dzetadx_s(i,N) = -1/H_s(i)*(1/(3*dx)*(-hB_s(i-1)-3*hB_s(i)+4*hB_s(i+1)) +...
%                    (N-1)*dzeta/(3*dx)*(-H_s(i-1)-3*H_s(i)+4*H_s(i+1)));
%    dzetadx_s(i,1) = -1/H_s(i)*(1/(3*dx)*(-hB_s(i-1)-3*hB_s(i)+4*hB_s(i+1)) +...
%                    (1-1)*dzeta/(3*dx)*(-H_s(i-1)-3*H_s(i)+4*H_s(i+1)));
%    elseif i == M_s
%    dzetadx_s(i,j) = -((hB_s(M_s)-hB_s(i-1))/(dx/2)+...
%                    (j-1)*dzeta*(H_s(M_s)-H_s(i-1))/(dx/2))/H_s(i);
%    else
%    dzetadx_s(i,j) = -((hB_s(i+1)-hB_s(i-1))/(2*dx)+...
%                    (j-1)*dzeta*(H_s(i+1)-H_s(i-1))/(2*dx))/H_s(i);
%    dzetadx_s(i,N) = -((hB_s(i+1)-hB_s(i-1))/(2*dx)+...
%                    (N-1)*dzeta*(H_s(i+1)-H_s(i-1))/(2*dx))/H_s(i);
%    dzetadx_s(i,1) = -((hB_s(i+1)-hB_s(i-1))/(2*dx)+...
%                    (1-1)*dzeta*(H_s(i+1)-H_s(i-1))/(2*dx))/H_s(i);
%    end
%    dzetadx_s(M_s,N) = -((hB_s(M_s)-hB_s(M_s-1))/(dx/2)+...
%                    (N-1)*dzeta*(H_s(M_s)-H_s(M_s-1))/(dx/2))/H_s(M_s);
%end
%end
%
%dzetadx = [dzetadx_s(1,:);(dzetadx_s(2:M_s-2,:)+dzetadx_s(3:M_s-1,:))/2;dzetadx_s(M_s,:)];

%save('ice_geo.mat','M','M_s','N','xi','dx','dzeta','dzetadx_s','zeta',...
%    'hB','hB_s','hS','H','H_s','dhSdx_s');
