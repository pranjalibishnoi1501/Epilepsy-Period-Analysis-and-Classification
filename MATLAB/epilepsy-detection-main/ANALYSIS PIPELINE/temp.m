N = 49; % Number of nodes
Layout = rand(N, 2); % Initialize layout matrix




% Plot the layout to visualize the node positions
L=20;
data = load('Seizure_Timesig_Sub_uni11.mat');

data = data.C;
num = size(data);
N = 500;
for count = 1:floor(N/L)
    disp(count)
    [adj_ring, adj_rand] = generateGraphs(49);
    W = adj_ring;
    top = adj_ring ~= 0;
    plotting(Layout, W, top, 1:49)
end

