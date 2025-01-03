% Initialize parameters
num_nodes = 33;
timesteps = 4267;
edgeProb = rand(num_nodes); % Random edge probabilities (adjust as needed)
seedSet = randperm(num_nodes, 5); % Random seed set with 5 nodes (adjust as needed)

% Initialize adjacency matrix over time
adjacency_matrices = zeros(num_nodes, num_nodes, timesteps);

% Run diffusion model for each timestep
for t = 1:timesteps
    % Initialize influence vector for current timestep
    influence = zeros(num_nodes, 1);
    
    % Initialize active set with seed nodes
    activeSet = seedSet;
    
    % Mark seed nodes as activated
    influence(seedSet) = 1;
    
    % Create adjacency matrix for current timestep
    A_timestep = zeros(num_nodes);
    
    while ~isempty(activeSet)
        % Pop a node from the active set
        u = activeSet(1);
        activeSet = activeSet(2:end);
        
        % Get neighbors of the current node
        neighbors = find(edgeProb(u, :) > 0);
        
        % Iterate over neighbors
        for v = neighbors
            % Check if neighbor is not activated and activate with edge probability
            if influence(v) == 0 && rand() <= edgeProb(u, v)
                influence(v) = 1;
                % Add neighbor to the active set
                activeSet = [activeSet, v];
                
                % Update adjacency matrix
                A_timestep(u, v) = 1;
                A_timestep(v, u) = 1; % Assuming undirected graph
            end
        end
    end
    
    % Store adjacency matrix for the current timestep
    adjacency_matrices(:,:,t) = A_timestep;
end

save('independent_cascade.mat',"adjacency_matrices");

% Display adjacency matrices (optional)
% To visualize adjacency matrix at a specific timestep (e.g., timestep 100)
% disp(adjacency_matrices(:,:,100));
