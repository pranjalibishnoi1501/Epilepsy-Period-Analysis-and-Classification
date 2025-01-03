% clc
% clear
% %graphs = load("theta_A_W_uni12_100ms.mat").W;
% graphs = load("theta_SeizurePartitioning.mat").W;
% % Number of nodes
% L = size(graphs, 1); % assuming adjacency_matrices is your 3D array
% % Number of time windows
% T = size(graphs, 3);
% 
% Parameters for the SIS model
beta = 0.3;  % infection rate
alpha = 0.1; % reactivation rate
gamma = 0.05; % rate of returning to susceptible state
% initial_infection_rate = 0.1; % initial fraction of infected nodes
% 
% % Initialize arrays to store results
% infected_fraction = zeros(L, T); % fraction of infected nodes at each time window
% 
% % Loop over time windows
% for tt = 1:T
%     % Get adjacency matrix for the current time window
%     A = graphs(:,:,tt);
% 
%     N = L;
%     num_nodes = L;
% 
%     % Initial conditions
%     S0 = 0.9 * ones(N, 1); % initial fraction of susceptible nodes
%     I0 = 0.1 * ones(N, 1); % initial fraction of infected nodes
%     Z0 = zeros(N, 1);      % initial fraction of zombified nodes
% 
%     % Combine initial conditions
%     y0 = [S0; I0; Z0];
% 
%     % Time vector
%     tspan = [0 100];
%     t = tspan(1):1:tspan(2);
% 
%     % Define ODE system
%     ode = @(t, y) [ -beta * y(1:N) .* (A * y(N+1:2*N)) + gamma * y(2*N+1:end);
%                      beta * y(1:N) .* (A * y(N+1:2*N)) - alpha * y(N+1:2*N);
%                      alpha * y(N+1:2*N) - gamma * y(2*N+1:end) ];
% 
%     % Solve ODEs
%     [t, y] = ode45(ode, tspan, y0);
% 
%     % y = y ./ sum(y,2);
% 
%     % Extract final states
%     final_states = y(end, :);
%     % normalizer = final_states(1:L) + final_states(L+1:2*L) + final_states(2*L+1:end);
%     normalizer = 1;
%     S_final = final_states(1:num_nodes) / normalizer;
%     I_final = final_states(num_nodes+1:2*num_nodes) / normalizer;
%     Z_final = final_states(2*num_nodes+1:end) / normalizer;
% 
%     infected_fraction(:,tt) = I_final;
% 
% 
% end
% 
% average_node_infection_rate = zeros(L,3);
% 
% % Calculate average infection fraction for each node
% avg_infected_fraction_first_250 = mean(infected_fraction(:, 1:250), 2);
% avg_infected_fraction_last_250 = mean(infected_fraction(:, end-249:end), 2);
% avg_infected_fraction_remaining = mean(infected_fraction(:, 251:end-250), 2);
% 
% average_node_infection_rate(:,1) = avg_infected_fraction_first_250;
% average_node_infection_rate(:,2) = avg_infected_fraction_remaining;
% average_node_infection_rate(:,3) = avg_infected_fraction_last_250;
% 
% graph_infection_rate = mean(average_node_infection_rate, 1);
% disp(graph_infection_rate)
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

% Initialize colors for scatter plots
colors = lines(8); % Change 7 to the number of datasets you have

% Loop over datasets
for file_index = 12:12 % Assuming datasets are named uni12, uni13, ..., uni18
    % Load data file
    datafile = sprintf("theta_A_W_uni%d_100ms.mat", file_index);
    data = load(datafile);
    graphs = data.W;
    
    % Number of nodes and time windows
    L = size(graphs, 1);
    T = size(graphs, 3);

    initial_infection_rate = 0.1; 

    % Initialize arrays to store results for this file
    infected_fraction = zeros(L, T); % fraction of infected nodes at each time window
    
    % Loop over time windows
    for tt = 1:T
        % Get adjacency matrix for the current time window
        A = graphs(:,:,tt);

        N = L;
        num_nodes = L;
    
        % Initial conditions
        S0 = 0.9 * ones(L, 1); % initial fraction of susceptible individuals
        I0 = 0.1 * ones(L, 1); % initial fraction of infected individuals
        Z0 = zeros(L, 1);      % initial fraction of zombified individuals
        y0 = [S0; I0; Z0];
    
        % Time vector
        tspan = [0 100];
        t = tspan(1):1:tspan(2);
    
        % ODE system
        ode = @(t, y) [ -beta * y(1:N) .* (A * y(N+1:2*N)) + gamma * y(2*N+1:end);
                      beta * y(1:N) .* (A * y(N+1:2*N)) - alpha * y(N+1:2*N);
                      alpha * y(N+1:2*N) - gamma * y(2*N+1:end) ];
    
        % Solve ODEs
        [~, y] = ode45(ode, tspan, y0);
    
        % Extract final states of each node
        final_states = y(end, :);
        normalizer = 1;
        I_final = final_states(L+1:2*L) / normalizer;
    
        infected_fraction(:,tt) = I_final;
    end
    
    % Calculate average infection fraction for each node
    avg_infected_fraction_first_250 = mean(infected_fraction(:, 1:250), 2);
    avg_infected_fraction_last_250 = mean(infected_fraction(:, end-249:end), 2);
    avg_infected_fraction_remaining = mean(infected_fraction(:, 251:end-250), 2);
    
    % Store results for this file
    average_node_infection_rate = zeros(L, 3);
    average_node_infection_rate(:,1) = avg_infected_fraction_first_250;
    average_node_infection_rate(:,2) = avg_infected_fraction_remaining;
    average_node_infection_rate(:,3) = avg_infected_fraction_last_250;
    
    % Plot scatter and fit line for this file
    x = average_node_infection_rate(:,1);
    y = average_node_infection_rate(:,2);
    coefficients = polyfit(x, y, 1);
    slope = coefficients(1);
    intercept = coefficients(2);
    color = colors(file_index - 10, :);
    scatter(x, y, 'o');
    hold on;
    x_fit = linspace(min(x), max(x), 100);
    y_fit = slope * x_fit + intercept;
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
end

% % Customize the plot
% xlabel('X Axis');
% ylabel('Y Axis');
% legend('Data 1', 'Linear Fit 1', 'Data 2', 'Linear Fit 2', 'Data 3', 'Linear Fit 3', 'Data 4', 'Linear Fit 4', 'Data 5', 'Linear Fit 5', 'Data 6', 'Linear Fit 6', 'Data 7', 'Linear Fit 7', 'Data 8', 'Linear Fit 8');
% grid on;
% 
% % Turn off the hold for future plots
% hold off;