% Performs static condensation and removes Dirichlet nodes from the global
% matrix solve Ka=F

function [K_condensed, F_condensed] = condensation(K, F, num_nodes, dirichlet_nodes)

K_condensed_rows = zeros(num_nodes - length(dirichlet_nodes(1, :)), num_nodes);
F_condensed = zeros(num_nodes - length(dirichlet_nodes(1, :)), 1);

K_row = 1;
i = 1;      % index for dirichlet_nodes
j = 1;      % index for condensed row

for K_row = 1:num_nodes
    if (find(dirichlet_nodes(1, :) == K_row))
        i = i + 1;
    else
        K_condensed_rows(j,:) = K(K_row,:);
        F_condensed(j) = F(j);
        j = j + 1;
    end
end

% perform static condensation to remove Dirichlet node columns from solve
K_condensed = zeros(num_nodes - length(dirichlet_nodes(1, :)), num_nodes - length(dirichlet_nodes(1, :)));

K_column = 1;
i = 1;          % index for dirichlet nodes
j = 1;          % index for condensed column
        
for K_column = 1:num_nodes
    if (find(dirichlet_nodes(1, :) == K_column))
        i = i + 1;
    else
        K_condensed(:,j) = K_condensed_rows(:, K_column);
        j = j + 1;
    end
end