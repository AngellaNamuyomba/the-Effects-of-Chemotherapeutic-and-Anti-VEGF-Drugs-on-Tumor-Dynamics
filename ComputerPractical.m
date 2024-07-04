%% COMPUTER PRACTICAL  

%% EXTENSION TO LOGISTIC GROWTH FUNCTION 
% Set the time of the simulation
t0 = 0; tf = 10;  %(day)

% Set the initial condition
% C, D, E, V, A
C0 = 49.0497;       %cancer cell volume
E0 = 1;             %endothelial cell density
V0 = 1.0;           %VEGF concentration  

% define drug doses
D_doses = [0.03, 0]; % doses for chemotherapy drug
A_doses = [0.05, 0]; % dosis for anti-VEGF drug
colors = lines(length(D_doses)*length(A_doses)); % generate unique colors

% initialize figure for plotting
figure;

%create subplots
subplot(2, 3, 1); hold on; title('Cancer cell volume(mm^3)'); xlabel('Time(day)');ylabel('Volume(mm^3');
subplot(2, 3, 2); hold on; title('Chemotherapeutic Drug Concentration (mg/ml)'); xlabel('time (day)'); ylabel('Concentration (mg/ml)');
subplot(2, 3, 3); hold on; title('Endothelial Cell Density (mm^3)'); xlabel('time (day)'); ylabel('Density (mm^3)');
subplot(2, 3, 4); hold on; title('Anti-VEGF Drug Concentration (mg/ml)'); xlabel('time (day)'); ylabel('Concentration (mg/ml)');
subplot(2, 3, 5); hold on; title('VEGF Concentration (mg/ml)'); xlabel('time (day)'); ylabel('Concentration (mg/ml)');

color_index = 1;

% loop over each combination of drug doses
for i = 1:length(D_doses)
    for j = 1:length(A_doses)
        D0 = D_doses(i); % current drug dose
        A0 = A_doses(j); % current drug dose
        x0 = [C0, D0, E0, A0, V0]; % initial condition


% Solve the ODE system
%opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[T,X] = ode15s(@Model3_reduced1,[t0 tf],x0,[]);

%plot each variable in it respective subplot
subplot(2, 3, 1);
        plot(T, X(:, 1), 'DisplayName', ['D0 = ', num2str(D0), ', A0 = ', num2str(A0)], 'Color', colors(color_index, :));
        
        subplot(2, 3, 2);
        plot(T, X(:, 2), 'DisplayName', ['D0 = ', num2str(D0), ', A0 = ', num2str(A0)], 'Color', colors(color_index, :));
        
        subplot(2, 3, 3);
        plot(T, X(:, 3), 'DisplayName', ['D0 = ', num2str(D0), ', A0 = ', num2str(A0)], 'Color', colors(color_index, :));
        
        subplot(2, 3, 4);
        plot(T, X(:, 4), 'DisplayName', ['D0 = ', num2str(D0), ', A0 = ', num2str(A0)], 'Color', colors(color_index, :));
        
        subplot(2, 3, 5);
        plot(T, X(:, 5), 'DisplayName', ['D0 = ', num2str(D0), ', A0 = ', num2str(A0)], 'Color', colors(color_index, :));
        
        color_index = color_index + 1;
    end
end

% Show legends
subplot(2, 3, 1); legend show;
subplot(2, 3, 2); legend show;
subplot(2, 3, 3); legend show;
subplot(2, 3, 4); legend show;
subplot(2, 3, 5); legend show;

subplot(2,3,5);
plot(T, X(:, 5));
xlabel('time (day)');
ylabel('anti-VEGF drug concentration B (mg/ml)');

