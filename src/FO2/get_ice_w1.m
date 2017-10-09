function [w_s, w] = get_ice_w1(u_s, u)

global rho g n kc Cp SPY ...
        M M_s N xi dx dzeta dzetadx dzetadx_s zeta hB hB_s hS H H_s dhSdx dhSdx_s dt ...
        de0 Sigma0


w = zeros(M,N);
w_s = zeros(M,N-1);
wc = zeros(1,N);
wc_sv = zeros(1,N+1);

for i = 1:M
    
    uc = u(i,:);
    uc_sv = [uc(1),(uc(1:N-1)+uc(2:N))/2,uc(N)];
        % vertical staggered horizontal v of column i
    for j = 1:N+1
        if j == 1 || j == 2 || j == N+1
            uc_bar_sv(i,1) = 0;
            uc_bar_sv(i,2) = uc_sv(2)/2;
            uc_bar_sv(i,N+1) = (sum(uc(2:N-1))+uc_sv(2)/4+(uc(N)+uc_sv(N))/4)/(N-1);
        else
            uc_bar_sv(i,j) = (sum(uc(2:j-1))+uc_sv(2)/4)/(j-3/2);
        end
    end
        % caculate vertical staggered u bar at j
end

for i = 1:M
    
    uc = u(i,:);

    if i == 1
        for j = 1:N+1
            if j == 1
                wc_sv(1,1) = 0;
            elseif j == N+1
                wc_sv(1,j) = -(j-2)*dzeta*(H(i+1)*uc_bar_sv(i+1,j)-H(i)*uc_bar_sv(i,j))/(dx)+...
                            (j-2)*dzeta*uc_sv(1,j)*(H(i+1)-H(i))/(dx) + uc_sv(1,j)*(hB(i+1)-hB(i))/(dx);
            else
                wc_sv(1,j) = -((j-1)-1/2)*dzeta*(H(i+1)*uc_bar_sv(i+1,j)-H(i)*uc_bar_sv(i,j))/(dx)+...
                            ((j-1)-1/2)*dzeta*uc_sv(1,j)*(H(i+1)-H(i))/(dx) + uc_sv(1,j)*(hB(i+1)-hB(i))/(dx);
            end
        end
    elseif i == M
        for j = 1:N+1
            if j == 1
                wc_sv(1,1) = 0;
            elseif j == N+1
                wc_sv(1,j) = -(j-2)*dzeta*(H(i)*uc_bar_sv(i,j)-H(i-1)*uc_bar_sv(i-1,j))/(dx)+...
                            (j-2)*dzeta*uc_sv(1,j)*(H(i)-H(i-1))/(dx) + uc_sv(1,j)*(hB(i)-hB(i-1))/(dx);
            else
                wc_sv(1,j) = -((j-1)-1/2)*dzeta*(H(i)*uc_bar_sv(i,j)-H(i-1)*uc_bar_sv(i-1,j))/(dx)+...
                            ((j-1)-1/2)*dzeta*uc_sv(1,j)*(H(i)-H(i-1))/(dx) + uc_sv(1,j)*(hB(i)-hB(i-1))/(dx);
            end
        end
    else
        for j = 1:N+1
            if j == 1
                wc_sv(1,1) = 0;
            elseif j == N+1
                wc_sv(1,j) = -(j-2)*dzeta*(H(i+1)*uc_bar_sv(i+1,j)-H(i-1)*uc_bar_sv(i-1,j))/(2*dx)+...
                            (j-2)*dzeta*uc_sv(1,j)*(H(i+1)-H(i-1))/(2*dx) + uc_sv(1,j)*(hB(i+1)-hB(i-1))/(2*dx);
            else
                wc_sv(1,j) = -((j-1)-1/2)*dzeta*(H(i+1)*uc_bar_sv(i+1,j)-H(i-1)*uc_bar_sv(i-1,j))/(2*dx)+...
                            ((j-1)-1/2)*dzeta*uc_sv(1,j)*(H(i+1)-H(i-1))/(2*dx) + uc_sv(1,j)*(hB(i+1)-hB(i-1))/(2*dx);
            end
        end
    end

    wc = [wc_sv(1), (wc_sv(2:N-1)+wc_sv(3:N))/2, wc_sv(N+1)];
    w(i,:) = wc;
    w_s(i,:) = wc_sv(2:end-1);

end

