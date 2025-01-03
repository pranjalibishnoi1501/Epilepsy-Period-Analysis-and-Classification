clc
clear

%datafile = "theta_A_W_uni20_100ms.mat";
%graphs = load(datafile).W;
%graphs = load("theta_A_W_uni20_100ms.mat").W;
% Number of nodes
%L = size(graphs, 1); % assuming adjacency_matrices is your 3D array
% Number of time windows
%T = size(graphs, 3);

% % Parameters for the SIS model
% beta = 1.0; % infection rate
% gamma = 0.1; % recovery rate
% initial_infection_rate = 0.1; % initial fraction of infected nodes

% Initialize arrays to store results
%infected_fraction = zeros(L, T); % fraction of infected nodes at each time window

% Loop over time windows
%for tt = 1:T
    % Get adjacency matrix for the current time window
    %A = graphs(:,:,tt);
    
    % Initial conditions
    %S0 = 0.9 * ones(L, 1); % initial fraction of susceptible individuals
    %I0 = 0.1 * ones(L, 1); % initial fraction of infected individuals
    %y0 = [S0; I0];
    
    % Time vector
    %tspan = [0 100];
    %t = tspan(1):1:tspan(2);
    
    % ODE system
    %ode = @(t, y) [gamma*y(L+1:end) - beta*(A*y(L+1:end)).*y(1:L);
                   %beta*(A*y(L+1:end)).*y(1:L) - gamma*y(L+1:end)];
    
    % Solve ODEs
    %[t, y] = ode45(ode, tspan, y0);

    %final_states = y(end, :);
    % normalizer = final_states(1:L) + final_states(L+1:end);
    %normalizer = 1;
    %S_final = final_states(1:L) / normalizer;
    %I_final = final_states(L+1:end) / normalizer;

    %infected_fraction(:,tt) = I_final;


%end

%average_node_infection_rate = zeros(L,3);

% Calculate average infection fraction for each node
% avg_infected_fraction_first_250 = mean(infected_fraction(:, 1:250), 2);
% avg_infected_fraction_last_250 = mean(infected_fraction(:, end-249:end), 2);
% avg_infected_fraction_remaining = mean(infected_fraction(:, 251:end-250), 2);
% 
% average_node_infection_rate(:,1) = avg_infected_fraction_first_250;
% average_node_infection_rate(:,2) = avg_infected_fraction_remaining;
% average_node_infection_rate(:,3) = avg_infected_fraction_last_250;
% 
% disp(average_node_infection_rate)

% % Define colors for scatter plots
% colors = {'b', 'g', 'r', 'c', 'm', 'y', 'k', [0.5 0.5 0.5]}; % You can add more colors as needed
% 
% % Define the number of datasets
% num_datasets = 8;
% 
% % Initialize arrays to store all data points
% all_x = [];
% all_y = [];
% 
% % Loop over each dataset
% for dataset = 1:num_datasets
%     % Load data
%     datafile = sprintf("theta_A_W_uni%d_100ms.mat", dataset + 10); % Adjust the index to match your filenames
%     data = load(datafile);
%     graphs = data.W;
%     L = size(graphs, 1);
%     T = size(graphs, 3);
% 
%     % Initialize arrays to store infected fractions
%     infected_fraction = zeros(L, T);
% 
%     % Calculate infected fraction for each time window
%     for tt = 1:T
%         A = graphs(:,:,tt);
%         S0 = 0.9 * ones(L, 1);
%         I0 = 0.1 * ones(L, 1);
%         y0 = [S0; I0];
%         tspan = [0 100];
%         t = tspan(1):1:tspan(2);
%         ode = @(t, y) [gamma*y(L+1:end) - beta*(A*y(L+1:end)).*y(1:L);
%                        beta*(A*y(L+1:end)).*y(1:L) - gamma*y(L+1:end)];
%         [t, y] = ode45(ode, tspan, y0);
%         final_states = y(end, :);
%         normalizer = 1;
%         S_final = final_states(1:L) / normalizer;
%         I_final = final_states(L+1:end) / normalizer;
%         infected_fraction(:,tt) = I_final;
%     end
% 
%     % Calculate average infected fraction
%     avg_infected_fraction_first_250 = mean(infected_fraction(:, 1:250), 2);
%     avg_infected_fraction_last_250 = mean(infected_fraction(:, end-249:end), 2);
%     avg_infected_fraction_remaining = mean(infected_fraction(:, 251:end-250), 2);
% 
%     % Store average infected fraction
%     average_node_infection_rate = [avg_infected_fraction_first_250, avg_infected_fraction_remaining, avg_infected_fraction_last_250];
% 
%     % Extract x and y values
%     x = average_node_infection_rate(:,1);
%     y = average_node_infection_rate(:,2);
% 
%     % Store x and y values
%     all_x = [all_x; x];
%     all_y = [all_y; y];
% 
%     % Plot data with different colors
%     scatter(x, y, [], colors{dataset}, 'o');
%     hold on;
% end
% 
% % Fit a linear model to all data points
% coefficients = polyfit(all_x, all_y, 1);
% slope = coefficients(1);
% intercept = coefficients(2);
% 
% % Plot the single linear fit line
% x_fit = linspace(min(all_x), max(all_x), 100);
% y_fit = slope * x_fit + intercept;
% plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
% 
% % Customize the plot
% % xlabel('Average Infected Fraction (First 250 Windows)');
% % ylabel('Average Infected Fraction (Remaining Windows)');
% title('Superimposed Scatter Plots with Linear Fit');
% 
% % Generate legend
% legend('Data 1', 'Data 2', 'Data 3', 'Data 4', 'Data 5', 'Data 6', 'Data 7', 'Data 8', 'Linear Fit');
% 
% grid on;
% hold off;

% Parameters for the SIS model
beta = 1.0; % infection rate
gamma = 0.1; % recovery rate
initial_infection_rate = 0.1; % initial fraction of infected nodes

% Load and plot the first dataset
datafile1 = "theta_A_W_uni12_100ms.mat";
data1 = load(datafile1);
graphs1 = data1.W;
L1 = size(graphs1, 1);
T1 = size(graphs1, 3);
infected_fraction1 = zeros(L1, T1);
% Get adjacency matrix for the current time window
tt = 150;
A = graphs1(:,:,tt);

% Initial conditions
S0 = 0.9 * ones(L1, 1); % initial fraction of susceptible individuals
I0 = 0.1 * ones(L1, 1); % initial fraction of infected individuals
y0 = [S0; I0];

% Time vector
tspan = [0 100];
t = tspan(1):1:tspan(2);

% ODE system
ode = @(t, y) [gamma*y(L1+1:end) - beta*(A*y(L1+1:end)).*y(1:L1);
               beta*(A*y(L1+1:end)).*y(1:L1) - gamma*y(L1+1:end)];

% Solve ODEs
[t, y] = ode45(ode, tspan, y0);

disp(size(y))

final_states = y(end, :);
% normalizer = final_states(1:L) + final_states(L+1:end);
normalizer = 1;
S_final = final_states(1:L1) / normalizer;
I_final = final_states(L1+1:end) / normalizer;

infected_fraction1(:,tt) = I_final;
% for tt = 1:T1
%     % Get adjacency matrix for the current time window
%     A = graphs1(:,:,tt);
% 
%     % Initial conditions
%     S0 = 0.9 * ones(L1, 1); % initial fraction of susceptible individuals
%     I0 = 0.1 * ones(L1, 1); % initial fraction of infected individuals
%     y0 = [S0; I0];
% 
%     % Time vector
%     tspan = [0 100];
%     t = tspan(1):1:tspan(2);
% 
%     % ODE system
%     ode = @(t, y) [gamma*y(L1+1:end) - beta*(A*y(L1+1:end)).*y(1:L1);
%                    beta*(A*y(L1+1:end)).*y(1:L1) - gamma*y(L1+1:end)];
% 
%     % Solve ODEs
%     [t, y] = ode45(ode, tspan, y0);
% 
%     final_states = y(end, :);
%     % normalizer = final_states(1:L) + final_states(L+1:end);
%     normalizer = 1;
%     S_final = final_states(1:L1) / normalizer;
%     I_final = final_states(L1+1:end) / normalizer;
% 
%     infected_fraction1(:,tt) = I_final;
% end
% 
% average_node_infection_rate1 = zeros(L1,3);
% 
% avg_infected_fraction_first1_250 = mean(infected_fraction1(:, 1:250), 2);
% avg_infected_fraction_last1_250 = mean(infected_fraction1(:, end-249:end), 2);
% avg_infected_fraction_remaining1 = mean(infected_fraction1(:, 251:end-250), 2);
% 
% average_node_infection_rate1(:,1) = avg_infected_fraction_first1_250;
% average_node_infection_rate1(:,2) = avg_infected_fraction_remaining1;
% average_node_infection_rate1(:,3) = avg_infected_fraction_last1_250;
% 
% disp(average_node_infection_rate1)
% x1 = average_node_infection_rate1(:,3);
% y1 = average_node_infection_rate1(:,2);
% coefficients1 = polyfit(x1, y1, 1);
% slope1 = coefficients1(1);
% intercept1 = coefficients1(2);
% scatter(x1, y1, 'o');
% hold on;
% x_fit1 = linspace(min(x1), max(x1), 100);
% y_fit1 = slope1 * x_fit1 + intercept1;
% plot(x_fit1, y_fit1, 'r-', 'LineWidth', 2);
% 
% % % Load and plot the second dataset
% % datafile2 = "theta_A_W_uni12_100ms.mat";
% % data2 = load(datafile2);
% % graphs2 = data2.W;
% % L2 = size(graphs2, 1);
% % T2 = size(graphs2, 3);
% % infected_fraction2 = zeros(L2, T2);
% % for tt = 1:T2
% %     % Get adjacency matrix for the current time window
% %     A = graphs2(:,:,tt);
% % 
% %     % Initial conditions
% %     S0 = 0.9 * ones(L2, 1); % initial fraction of susceptible individuals
% %     I0 = 0.1 * ones(L2, 1); % initial fraction of infected individuals
% %     y0 = [S0; I0];
% % 
% %     % Time vector
% %     tspan = [0 100];
% %     t = tspan(1):1:tspan(2);
% % 
% %     % ODE system
% %     ode = @(t, y) [gamma*y(L2+1:end) - beta*(A*y(L2+1:end)).*y(1:L2);
% %                    beta*(A*y(L2+1:end)).*y(1:L2) - gamma*y(L2+1:end)];
% % 
% %     % Solve ODEs
% %     [t, y] = ode45(ode, tspan, y0);
% % 
% %     final_states = y(end, :);
% %     % normalizer = final_states(1:L) + final_states(L+1:end);
% %     normalizer = 1;
% %     S_final = final_states(1:L2) / normalizer;
% %     I_final = final_states(L2+1:end) / normalizer;
% % 
% %     infected_fraction2(:,tt) = I_final;
% % end
% % 
% % average_node_infection_rate2 = zeros(L2,3);
% % 
% % avg_infected_fraction_first2_250 = mean(infected_fraction2(:, 1:250), 2);
% % avg_infected_fraction_last2_250 = mean(infected_fraction2(:, end-249:end), 2);
% % avg_infected_fraction_remaining2 = mean(infected_fraction2(:, 251:end-250), 2);
% % 
% % average_node_infection_rate2(:,1) = avg_infected_fraction_first2_250;
% % average_node_infection_rate2(:,2) = avg_infected_fraction_remaining2;
% % average_node_infection_rate2(:,3) = avg_infected_fraction_last2_250;
% % 
% % disp(average_node_infection_rate2)
% % x2 = average_node_infection_rate2(:,3);
% % y2 = average_node_infection_rate2(:,2);
% % coefficients2 = polyfit(x2, y2, 1);
% % slope2 = coefficients2(1);
% % intercept2 = coefficients2(2);
% % scatter(x2, y2, 'o', 'MarkerEdgeColor', 'g');
% % x_fit2 = linspace(min(x2), max(x2), 100);
% % y_fit2 = slope2 * x_fit2 + intercept2;
% % plot(x_fit2, y_fit2, 'r-', 'LineWidth', 2);
% % 
% % % Load and plot the third dataset
% % datafile3 = "theta_A_W_uni13_100ms.mat";
% % data3 = load(datafile3);
% % graphs3 = data3.W;
% % L3 = size(graphs3, 1);
% % T3 = size(graphs3, 3);
% % infected_fraction3 = zeros(L3, T3);
% % for tt = 1:T3
% %     % Get adjacency matrix for the current time window
% %     A = graphs3(:,:,tt);
% % 
% %     % Initial conditions
% %     S0 = 0.9 * ones(L3, 1); % initial fraction of susceptible individuals
% %     I0 = 0.1 * ones(L3, 1); % initial fraction of infected individuals
% %     y0 = [S0; I0];
% % 
% %     % Time vector
% %     tspan = [0 100];
% %     t = tspan(1):1:tspan(2);
% % 
% %     % ODE system
% %     ode = @(t, y) [gamma*y(L3+1:end) - beta*(A*y(L3+1:end)).*y(1:L3);
% %                    beta*(A*y(L3+1:end)).*y(1:L3) - gamma*y(L3+1:end)];
% % 
% %     % Solve ODEs
% %     [t, y] = ode45(ode, tspan, y0);
% % 
% %     final_states = y(end, :);
% %     % normalizer = final_states(1:L) + final_states(L+1:end);
% %     normalizer = 1;
% %     S_final = final_states(1:L3) / normalizer;
% %     I_final = final_states(L3+1:end) / normalizer;
% % 
% %     infected_fraction3(:,tt) = I_final;
% % end
% % 
% % average_node_infection_rate3 = zeros(L3,3);
% % 
% % avg_infected_fraction_first3_250 = mean(infected_fraction3(:, 1:250), 2);
% % avg_infected_fraction_last3_250 = mean(infected_fraction3(:, end-249:end), 2);
% % avg_infected_fraction_remaining3 = mean(infected_fraction3(:, 251:end-250), 2);
% % 
% % average_node_infection_rate3(:,1) = avg_infected_fraction_first3_250;
% % average_node_infection_rate3(:,2) = avg_infected_fraction_remaining3;
% % average_node_infection_rate3(:,3) = avg_infected_fraction_last3_250;
% % 
% % disp(average_node_infection_rate3)
% % x3 = average_node_infection_rate3(:,3);
% % y3 = average_node_infection_rate3(:,2);
% % coefficients3 = polyfit(x3, y3, 1);
% % slope3 = coefficients3(1);
% % intercept3 = coefficients3(2);
% % scatter(x3, y3, 'o', 'MarkerEdgeColor', 'r');
% % hold on;
% % x_fit3 = linspace(min(x3), max(x3), 100);
% % y_fit3 = slope3 * x_fit3 + intercept3;
% % plot(x_fit3, y_fit3, 'r-', 'LineWidth', 2);
% % 
% % % Load and plot the fourth dataset
% % datafile4 = "theta_A_W_uni14_100ms.mat";
% % data4 = load(datafile4);
% % graphs4 = data4.W;
% % L4 = size(graphs4, 1);
% % T4 = size(graphs4, 3);
% % infected_fraction4 = zeros(L4, T4);
% % for tt = 1:T4
% %     % Get adjacency matrix for the current time window
% %     A = graphs4(:,:,tt);
% % 
% %     % Initial conditions
% %     S0 = 0.9 * ones(L4, 1); % initial fraction of susceptible individuals
% %     I0 = 0.1 * ones(L4, 1); % initial fraction of infected individuals
% %     y0 = [S0; I0];
% % 
% %     % Time vector
% %     tspan = [0 100];
% %     t = tspan(1):1:tspan(2);
% % 
% %     % ODE system
% %     ode = @(t, y) [gamma*y(L4+1:end) - beta*(A*y(L4+1:end)).*y(1:L4);
% %                    beta*(A*y(L4+1:end)).*y(1:L4) - gamma*y(L4+1:end)];
% % 
% %     % Solve ODEs
% %     [t, y] = ode45(ode, tspan, y0);
% % 
% %     final_states = y(end, :);
% %     % normalizer = final_states(1:L) + final_states(L+1:end);
% %     normalizer = 1;
% %     S_final = final_states(1:L4) / normalizer;
% %     I_final = final_states(L4+1:end) / normalizer;
% % 
% %     infected_fraction4(:,tt) = I_final;
% % end
% % 
% % average_node_infection_rate4 = zeros(L4,3);
% % 
% % avg_infected_fraction_first4_250 = mean(infected_fraction4(:, 1:250), 2);
% % avg_infected_fraction_last4_250 = mean(infected_fraction4(:, end-249:end), 2);
% % avg_infected_fraction_remaining4 = mean(infected_fraction4(:, 251:end-250), 2);
% % 
% % average_node_infection_rate4(:,1) = avg_infected_fraction_first4_250;
% % average_node_infection_rate4(:,2) = avg_infected_fraction_remaining4;
% % average_node_infection_rate4(:,3) = avg_infected_fraction_last4_250;
% % 
% % disp(average_node_infection_rate4)
% % x4 = average_node_infection_rate4(:,3);
% % y4 = average_node_infection_rate4(:,2);
% % coefficients4 = polyfit(x4, y4, 1);
% % slope4 = coefficients4(1);
% % intercept4 = coefficients4(2);
% % scatter(x4, y4, 'o', 'MarkerEdgeColor', 'm');
% % hold on;
% % x_fit4 = linspace(min(x4), max(x4), 100);
% % y_fit4 = slope4 * x_fit4 + intercept4;
% % plot(x_fit4, y_fit4, 'r-', 'LineWidth', 2);
% % 
% % % Load and plot the fifth dataset
% % datafile5 = "theta_A_W_uni15_100ms.mat";
% % data5 = load(datafile5);
% % graphs5 = data5.W;
% % L5 = size(graphs5, 1);
% % T5 = size(graphs5, 3);
% % infected_fraction5 = zeros(L5, T5);
% % for tt = 1:T5
% %     % Get adjacency matrix for the current time window
% %     A = graphs5(:,:,tt);
% % 
% %     % Initial conditions
% %     S0 = 0.9 * ones(L5, 1); % initial fraction of susceptible individuals
% %     I0 = 0.1 * ones(L5, 1); % initial fraction of infected individuals
% %     y0 = [S0; I0];
% % 
% %     % Time vector
% %     tspan = [0 100];
% %     t = tspan(1):1:tspan(2);
% % 
% %     % ODE system
% %     ode = @(t, y) [gamma*y(L5+1:end) - beta*(A*y(L5+1:end)).*y(1:L5);
% %                    beta*(A*y(L5+1:end)).*y(1:L5) - gamma*y(L5+1:end)];
% % 
% %     % Solve ODEs
% %     [t, y] = ode45(ode, tspan, y0);
% % 
% %     final_states = y(end, :);
% %     % normalizer = final_states(1:L) + final_states(L+1:end);
% %     normalizer = 1;
% %     S_final = final_states(1:L5) / normalizer;
% %     I_final = final_states(L5+1:end) / normalizer;
% % 
% %     infected_fraction5(:,tt) = I_final;
% % end
% % 
% % average_node_infection_rate5 = zeros(L5,3);
% % 
% % avg_infected_fraction_first5_250 = mean(infected_fraction5(:, 1:250), 2);
% % avg_infected_fraction_last5_250 = mean(infected_fraction5(:, end-249:end), 2);
% % avg_infected_fraction_remaining5 = mean(infected_fraction5(:, 251:end-250), 2);
% % 
% % average_node_infection_rate5(:,1) = avg_infected_fraction_first5_250;
% % average_node_infection_rate5(:,2) = avg_infected_fraction_remaining5;
% % average_node_infection_rate5(:,3) = avg_infected_fraction_last5_250;
% % 
% % disp(average_node_infection_rate5)
% % x5 = average_node_infection_rate5(:,3);
% % y5 = average_node_infection_rate5(:,2);
% % coefficients5 = polyfit(x5, y5, 1);
% % slope5 = coefficients5(1);
% % intercept5 = coefficients5(2);
% % scatter(x5, y5, 'o', 'MarkerEdgeColor', 'r');
% % hold on;
% % x_fit5 = linspace(min(x5), max(x5), 100);
% % y_fit5 = slope5 * x_fit5 + intercept5;
% % plot(x_fit5, y_fit5, 'r-', 'LineWidth', 2);
% % 
% % % Load and plot the sixth dataset
% % datafile6 = "theta_A_W_uni16_100ms.mat";
% % data6 = load(datafile6);
% % graphs6 = data6.W;
% % L6 = size(graphs6, 1);
% % T6 = size(graphs6, 3);
% % infected_fraction6 = zeros(L6, T6);
% % for tt = 1:T6
% %     % Get adjacency matrix for the current time window
% %     A = graphs6(:,:,tt);
% % 
% %     % Initial conditions
% %     S0 = 0.9 * ones(L6, 1); % initial fraction of susceptible individuals
% %     I0 = 0.1 * ones(L6, 1); % initial fraction of infected individuals
% %     y0 = [S0; I0];
% % 
% %     % Time vector
% %     tspan = [0 100];
% %     t = tspan(1):1:tspan(2);
% % 
% %     % ODE system
% %     ode = @(t, y) [gamma*y(L6+1:end) - beta*(A*y(L6+1:end)).*y(1:L6);
% %                    beta*(A*y(L6+1:end)).*y(1:L6) - gamma*y(L6+1:end)];
% % 
% %     % Solve ODEs
% %     [t, y] = ode45(ode, tspan, y0);
% % 
% %     final_states = y(end, :);
% %     % normalizer = final_states(1:L) + final_states(L+1:end);
% %     normalizer = 1;
% %     S_final = final_states(1:L6) / normalizer;
% %     I_final = final_states(L6+1:end) / normalizer;
% % 
% %     infected_fraction6(:,tt) = I_final;
% % end
% % 
% % average_node_infection_rate6 = zeros(L6,3);
% % 
% % avg_infected_fraction_first6_250 = mean(infected_fraction6(:, 1:250), 2);
% % avg_infected_fraction_last6_250 = mean(infected_fraction6(:, end-249:end), 2);
% % avg_infected_fraction_remaining6 = mean(infected_fraction6(:, 251:end-250), 2);
% % 
% % average_node_infection_rate6(:,1) = avg_infected_fraction_first6_250;
% % average_node_infection_rate6(:,2) = avg_infected_fraction_remaining6;
% % average_node_infection_rate6(:,3) = avg_infected_fraction_last6_250;
% % 
% % disp(average_node_infection_rate6)
% % x6 = average_node_infection_rate6(:,3);
% % y6 = average_node_infection_rate6(:,2);
% % coefficients6 = polyfit(x6, y6, 1);
% % slope6 = coefficients6(1);
% % intercept6 = coefficients6(2);
% % scatter(x6, y6, 'o', 'MarkerEdgeColor', 'c');
% % hold on;
% % x_fit6 = linspace(min(x6), max(x6), 100);
% % y_fit6 = slope6 * x_fit6 + intercept6;
% % plot(x_fit6, y_fit6, 'r-', 'LineWidth', 2);
% % 
% % % Load and plot the seventh dataset
% % datafile7 = "theta_A_W_uni17_100ms.mat";
% % data7 = load(datafile7);
% % graphs7 = data7.W;
% % L7 = size(graphs7, 1);
% % T7 = size(graphs7, 3);
% % infected_fraction7 = zeros(L7, T7);
% % for tt = 1:T7
% %     % Get adjacency matrix for the current time window
% %     A = graphs7(:,:,tt);
% % 
% %     % Initial conditions
% %     S0 = 0.9 * ones(L7, 1); % initial fraction of susceptible individuals
% %     I0 = 0.1 * ones(L7, 1); % initial fraction of infected individuals
% %     y0 = [S0; I0];
% % 
% %     % Time vector
% %     tspan = [0 100];
% %     t = tspan(1):1:tspan(2);
% % 
% %     % ODE system
% %     ode = @(t, y) [gamma*y(L7+1:end) - beta*(A*y(L7+1:end)).*y(1:L7);
% %                    beta*(A*y(L7+1:end)).*y(1:L7) - gamma*y(L7+1:end)];
% % 
% %     % Solve ODEs
% %     [t, y] = ode45(ode, tspan, y0);
% % 
% %     final_states = y(end, :);
% %     % normalizer = final_states(1:L) + final_states(L+1:end);
% %     normalizer = 1;
% %     S_final = final_states(1:L7) / normalizer;
% %     I_final = final_states(L7+1:end) / normalizer;
% % 
% %     infected_fraction7(:,tt) = I_final;
% % end
% % 
% % average_node_infection_rate7 = zeros(L7,3);
% % 
% % avg_infected_fraction_first7_250 = mean(infected_fraction7(:, 1:250), 2);
% % avg_infected_fraction_last7_250 = mean(infected_fraction7(:, end-249:end), 2);
% % avg_infected_fraction_remaining7 = mean(infected_fraction7(:, 251:end-250), 2);
% % 
% % average_node_infection_rate7(:,1) = avg_infected_fraction_first7_250;
% % average_node_infection_rate7(:,2) = avg_infected_fraction_remaining7;
% % average_node_infection_rate7(:,3) = avg_infected_fraction_last7_250;
% % 
% % disp(average_node_infection_rate7)
% % x7 = average_node_infection_rate7(:,3);
% % y7 = average_node_infection_rate7(:,2);
% % coefficients7 = polyfit(x7, y7, 1);
% % slope7 = coefficients7(1);
% % intercept7 = coefficients7(2);
% % scatter(x7, y7, 'o', 'MarkerEdgeColor', 'm');
% % hold on;
% % x_fit7 = linspace(min(x7), max(x7), 100);
% % y_fit7 = slope7 * x_fit7 + intercept7;
% % plot(x_fit7, y_fit7, 'r-', 'LineWidth', 2);
% % 
% % % Load and plot the eighth dataset
% % datafile8 = "theta_A_W_uni18_100ms.mat";
% % data8 = load(datafile8);
% % graphs8 = data8.W;
% % L8 = size(graphs8, 1);
% % T8 = size(graphs8, 3);
% % infected_fraction8 = zeros(L8, T8);
% % for tt = 1:T8
% %     % Get adjacency matrix for the current time window
% %     A = graphs8(:,:,tt);
% % 
% %     % Initial conditions
% %     S0 = 0.9 * ones(L8, 1); % initial fraction of susceptible individuals
% %     I0 = 0.1 * ones(L8, 1); % initial fraction of infected individuals
% %     y0 = [S0; I0];
% % 
% %     % Time vector
% %     tspan = [0 100];
% %     t = tspan(1):1:tspan(2);
% % 
% %     % ODE system
% %     ode = @(t, y) [gamma*y(L8+1:end) - beta*(A*y(L8+1:end)).*y(1:L8);
% %                    beta*(A*y(L8+1:end)).*y(1:L8) - gamma*y(L8+1:end)];
% % 
% %     % Solve ODEs
% %     [t, y] = ode45(ode, tspan, y0);
% % 
% %     final_states = y(end, :);
% %     % normalizer = final_states(1:L) + final_states(L+1:end);
% %     normalizer = 1;
% %     S_final = final_states(1:L8) / normalizer;
% %     I_final = final_states(L8+1:end) / normalizer;
% % 
% %     infected_fraction8(:,tt) = I_final;
% % end
% % 
% % average_node_infection_rate8 = zeros(L8,3);
% % 
% % avg_infected_fraction_first8_250 = mean(infected_fraction8(:, 1:250), 2);
% % avg_infected_fraction_last8_250 = mean(infected_fraction8(:, end-249:end), 2);
% % avg_infected_fraction_remaining8 = mean(infected_fraction8(:, 251:end-250), 2);
% % 
% % average_node_infection_rate8(:,1) = avg_infected_fraction_first8_250;
% % average_node_infection_rate8(:,2) = avg_infected_fraction_remaining8;
% % average_node_infection_rate8(:,3) = avg_infected_fraction_last8_250;
% % 
% % disp(average_node_infection_rate8)
% % x8 = average_node_infection_rate8(:,3);
% % y8 = average_node_infection_rate8(:,2);
% % coefficients8 = polyfit(x8, y8, 1);
% % slope8 = coefficients8(1);
% % intercept8 = coefficients8(2);
% % scatter(x8, y8, 'o', 'MarkerEdgeColor', 'y');
% % hold on;
% % x_fit8 = linspace(min(x8), max(x8), 100);
% % y_fit8 = slope8 * x_fit8 + intercept8;
% % plot(x_fit8, y_fit8, 'r-', 'LineWidth', 2);
% % 
% % % Customize the plot (labels, legend, etc.) if needed
% % xlabel('X-axis');
% % ylabel('Y-axis');
% % title('Superimposed Scatter Plots with Linear Fits');
% % legend('Data 1', 'Linear Fit 1', 'Data 2', 'Linear Fit 2', 'Data 3', 'Linear Fit 3', 'Data 4', 'Linear Fit 4', 'Data 5', 'Linear Fit 5', 'Data 6', 'Linear Fit 6', 'Data 7', 'Linear Fit 7', 'Data 8', 'Linear Fit 8');
% % grid on;
% % 
% % % Turn off the hold for future plots
% % hold off;
% % 
% % 
% % coefficients = polyfit(average_node_infection_rate(:,1), average_node_infection_rate(:,2), 1);
% % 
% % slope = coefficients(1);
% % intercept = coefficients(2);
% % 
% % scatter(average_node_infection_rate(:,1), average_node_infection_rate(:,2), 'o');
% % hold on;
% % 
% % x_fit = linspace(min(average_node_infection_rate(:,1)), max(average_node_infection_rate(:,1)), 100);
% % y_fit = slope * x_fit + intercept;
% % plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
% % 
% % grid on;
% % 
% % hold off;
% 
% graph_infection_rate = mean(average_node_infection_rate, 1);
% disp(graph_infection_rate)
% 
% save("sis_"+datafile,"average_node_infection_rate");
% 
% % Calculate correlation between columns of average_node_infection_rate
% correlation_matrix = corrcoef(average_node_infection_rate);
% 
% % Display correlation matrix
% disp('Correlation matrix between columns of average_node_infection_rate:');
% disp(correlation_matrix);
% 
% % Number of columns (time windows)
% num_columns = size(average_node_infection_rate, 2);
% 
% % Initialize matrices to store similarities
% cosine_similarity_matrix = zeros(num_columns, num_columns);
% euclidean_distance_matrix = zeros(num_columns, num_columns);
% 
% % Calculate cosine similarity and Euclidean distance between columns
% for i = 1:num_columns
%     for j = 1:num_columns
%         vector1 = average_node_infection_rate(:, i);
%         vector2 = average_node_infection_rate(:, j);
% 
%         % Cosine Similarity
%         cosine_similarity = dot(vector1, vector2) / (norm(vector1) * norm(vector2));
%         cosine_similarity_matrix(i, j) = cosine_similarity;
% 
%         % Euclidean Distance
%         euclidean_distance = norm(vector1 - vector2);
%         euclidean_distance_matrix(i, j) = euclidean_distance;
%     end
% end
% 
% % Display matrices
% disp('Cosine Similarity Matrix:');
% disp(cosine_similarity_matrix);
% disp('Euclidean Distance Matrix:');
% disp(euclidean_distance_matrix);
% 
% 
% % % Display average infection fraction for each node
% % for node = 1:L
% %     fprintf('Node %d: Average Infected Fraction (First 250 windows) = %.4f\n', node, avg_infected_fraction_first_250(node));
% %     fprintf('Node %d: Average Infected Fraction (Last 250 windows) = %.4f\n', node, avg_infected_fraction_last_250(node));
% %     fprintf('Node %d: Average Infected Fraction (Remaining windows) = %.4f\n', node, avg_infected_fraction_remaining(node));
% % end
% 
% % % Plot results
% % figure;
% % plot(1:T, infected_fraction, '-o');
% % xlabel('Time Window');
% % ylabel('Fraction of Infected Nodes');
% % title('SIS Model Simulation');
% % grid on;
