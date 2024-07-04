function dYdt = Model3_reduced(t,Y)
% ODE cancer model, including cancer cells (C), chemotherapy drug (e.g. docetaxel) (D)
%
% input: 
% - t:  time
% - Y:  vector of state variables (C, D)
% - dYdt:   column vector of derivatives

%% Parameters
c_M = 10^4; %  carrying capacity tumor cells mm^3
d_c = 0.1030;  % natural death rate tumor cells day-1
b = 0.1685; % death rate tumor cells by chemotherapy ml/mg.day
lambda_ce =0.4579; % proliferation rate tumor cell 1/day
d_d = 0.1825;  % clearance chemotherapeutic drug  1/day
b_k = 5*10^-5; % clearance chemotherapeutic drug binding cancer cells 1/day.mm^3
lambda_e = 0.03; % endothelial proliferation 1/day
d_e = 0.05; % natural death rate endothelial cells 1/day
k = 1; % rate VEGF induced proliferation 1/day
d_b = 0.05; % clearance anti-VEGF drug 1/day
u_bk = 0.15; % clearance anti-VEGF bound to VEGF ml/day.mg
lambda_v0 = 0.001; % baseline VEGF secretion by cancer cells mg/ml.day.mm
d_v = 0.1; % decay rate VEGF 1/day
u_b = 10; % removal VEGF bound antiVEGF ml/day.mg
e_M = 10;  % carrying capacity endothelial cells mm^3
k_m = 0.06 % micaelis-menten constant

%% Variables 
C = Y(1); % cancer cells 
D = Y(2); % chemotherapy drug concentration (e.g. docetaxel) 
E = Y(3); % endothelial density
A = Y(4); % anti-VEGF concentration
V = Y(5); % VEGF concentration
%5 variables -> 5 ODEs 

%% ODEs 

dCdt = lambda_ce*C*(1-C/c_M)-d_c*C-(b*C*D)/((k_m/E)+D);

dDdt = -d_d*D-b_k*C*D;

dEdt = lambda_e*E*(1-E/e_M)-d_e*E+k*V*E*(1-E/e_M);

dAdt = -d_b*A-u_bk*A*V;

dVdt = lambda_v0*C-d_v*V-u_b*V*A;

%create output column vector dYdt
dYdt = [dCdt;  dDdt; dEdt; dAdt; dVdt]; %5 variables ->5 ODEs

end

