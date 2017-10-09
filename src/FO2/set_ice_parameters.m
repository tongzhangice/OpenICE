function set_ice_parameters()

global rho g n kc Cp SPY ...
        M M_s N xi dx dzeta dzetadx dzetadx_s zeta hB hB_s hS H H_s dhSdx dhSdx_s dt ...
        de0 Sigma0

rho = 910;
    % density of ice [kg/m3]
g = 9.8;    
    % gravitational acceleration [m/s2]
n = 3;    
    % flow law exponent [ ]
kc = 2.1; 
    % conductivity of ice [W/m/K]
Cp = 2009;
    % specific heat capacity of ice [J/Kg/K]
de0 = 1e-60; Sigma0 = 1e-10;
    % small value in case of singularity
SPY = 31536000;

