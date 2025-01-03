function [adj_ring, adj_rand, largest_component_nodes] = generateGraphsAndLargestComponent(N)

%% Ring graph

adj_ring = zeros(N,N);
for i = 1:(N-1)
    adj_ring(i, i+1) = 1;
    adj_ring(i+1, i) = 1;
end
adj_ring(1,N) = 1;
adj_ring(N,1) = 1;

%% Random graph

adj_rand = adjProb(N,0.5);

%% Finding the largest connected component in the random graph
G_rand = graph(adj_rand);
bins = conncomp(G_rand);
component_sizes = histcounts(bins, 'BinMethod', 'integers');
[~, idx] = max(component_sizes);
largest_component_nodes = find(bins == idx);

end

function adj = adjProb(N, p)
    adj = rand(N) < p;
    adj = triu(adj,1);
    adj = adj + adj';
end
