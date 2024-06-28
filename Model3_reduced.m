function dYdt = Model3_reduced(t,Y)
% ODE cancer model, including cancer cells (C), chemotherapy drug (e.g. docetaxel) (D)
%
% input: 
% - t:  time
% - Y:  vector of state variables (C, D)
% - dYdt:   column vector of derivatives

%% Parameters
c_M = 10^4; % mm^3
d_c = 0.1030;  % day-1
b = 0.1685; % ml/mg.day
lambda_ce =0.4579; %1/day
d_d = 0.1825; %1/day
b_k = 5*10^-5; % % 1/day.mm^3
lambda_e = 0.03; 
d_e = 0.05; 
k = 1; 
d_b = 0.05; 
u_bk = 0.15; 
lambda_v0 = 0.001; 
E_thresh = 4.5; 
d_v = 0.1; 
u_b = 10; 
e_M = 10;

%% Variables 
C = Y(1); %cancer cells 
D = Y(2); %chemotherapy drug concentration (e.g. docetaxel) 
E = Y(3);
A = Y(4);
V = Y(5);
%2 variables -> 2 ODEs 

%% ODEs 

dCdt = lambda_ce*C*(1-C/c_M)-d_c*C-(b*C*D)/((k_m/E)+D);

dDdt = -d_d*D-b_k*C*D;

dEdt = lambda_e*E*(1-E/e_M)-d_e*E+k*V*E*(1-E/e_M);

dAdt = -d_b*A-u_bk*A*V;

dVdt = lambda_v0*C-d_v*V-u_b*V*A;

%create output column vector dYdt
dYdt = [dCdt;  dDdt; dEdt; dAdt; dVdt]; %2 variables -> 2 ODEs

end

