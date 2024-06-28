%% COMPUTER PRACTICAL  
%% EXTENSION TO LOGISTIC GROWTH FUNCTION 
% Set the time of the simulation
t0 = 0; tf = 19;  %(day)

% Set the initial condition
% C, D
x0 = [49.0497; 0.171; 0.5; 0.003; 0.5];  

% Solve the ODE system
%opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[T,X] = ode15s'Model3_reduced[1],[t0 tf];x0;[]';

%figure 2 plots
figure;
subplot(2, 3, 1);
plot(T, X(:, 1));
xlabel('time (day)');
ylabel('cancer cell volume (mm^3)');

subplot(2,3,2);
plot(T, X(:, 2));
xlabel('time (day)');
ylabel('chemotherapeutic drug concentration D (mg/ml)');

subplot(2,3,3);
plot(T, X(:, 3));
xlabel('time (day)');
ylabel('endothelial cell density E (mm^3)');

subplot(2,3,4);
plot(T, X(:, 4));
xlabel('time (day)');
ylabel('VEGF concentration V (mg/ml)');

subplot(2,3,5);
plot(T, X(:, 5));
xlabel('time (day)');
ylabel('anti-VEGF drug concentration B (mg/ml)');

