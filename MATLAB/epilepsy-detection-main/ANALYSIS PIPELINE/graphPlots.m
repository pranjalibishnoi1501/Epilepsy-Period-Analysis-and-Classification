clc
clear
%% Read adj
graphs = load("theta_A_W_uni11_100ms.mat").A;
graphs_w = load("theta_A_W_uni11_100ms.mat").W;
graphs_top = load("theta_A_W_uni11_100ms.mat").Top;

%graphs = load("theta_SeizurePartitioning.mat").A;
%graphs_w = load("theta_SeizurePartitioning.mat").W;
%graphs_top = load("theta_SeizurePartitioning.mat").Top;

N = size(graphs);
L = N(1);
N = N(3);
gc = zeros(N);
E = zeros(L,N);
C = zeros(L,N);
for i = 1:N
    adj = graphs(:,:,i);
    adj_w = graphs_w(:,:,i);
    [GC,gc_nodes]=giant_component(adj);
    %comp_mat = find_conn_comp(adj);
    %disp(size(comp_mat))
    gc(i) = length(gc_nodes);
    %gc(i) = length(comp_mat);
    %%%% Eigen centrality
    eigC = eigencentrality(adj_w);

    %%% Closeness Centrality
    clC = closeness(adj_w);

    E(:,i) = eigC;
    C(:,i) = clC;
end
gc1 = mean(gc(1:250));
gc2 = mean(gc(251:N-250));
gc3 = mean(gc(N-249:N));


E_compare = zeros(L,3);
C_compare = zeros(L,3);

e1 = mean(E(:,1:250),2);
e2 = mean(E(:,251:N-250),2);
e3 = mean(E(:,N-249:N),2);

E_compare(:,1) = e1;
E_compare(:,2) = e2;
E_compare(:,3) = e3;

c1 = mean(C(:,1:250),2);
c2 = mean(C(:,251:N-250),2);
c3 = mean(C(:,N-249:N),2);

C_compare(:,1) = c1;
C_compare(:,2) = c2;
C_compare(:,3) = c3;

disp(gc1)
disp(gc2)
disp(gc3)

disp(E_compare)
disp(C_compare)

%Layout = rand(L,2);
%plotting(Layout, adj1, adj1, 1:L)
%plotting(Layout, adj2, adj2_top, 1:L)


%% Run loop
%Layout = rand(L,2);
% N = 5000;
%for i = 1:N
%    adj = graphs(:,:,i);
%    plotting(Layout, adj, adj, 1:L)
%end