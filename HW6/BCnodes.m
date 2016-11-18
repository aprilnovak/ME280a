% Script to return the node numbers associated with different types of
% boundary conditions

function [dirichlet_nodes, neumann_nodes, a_k] = BCnodes(left, right, left_value, right_value, num_nodes, Nr, No)

% arrays that hold the nodes in the first row and the values in each column
dirichlet_nodes = [];
neumann_nodes = [];

% assign the nodes to either dirichlet or neumann BCs
i = 1;
switch left
    case 'Dirichlet'
        dirichlet_nodes(1, :) = 1:1:(Nr + 1);
        dirichlet_nodes(2, :) = left_value .* ones(1, Nr + 1);
    case 'Neumann'
        neumann_nodes(1, :) = 1:1:(Nr + 1);
        neumann_nodes(2,:) = left_value .* ones(1, Nr + 1);
    otherwise
        disp('You entered an incorrect field for the type of BC on the boundaries.');
end

i = i + 1;
switch right
    case 'Dirichlet'
        dirichlet_nodes(1, :) = (num_nodes - Nr):1:num_nodes;
        dirichlet_nodes(2, :) = right_value .* ones(1, Nr + 1);
    case 'Neumann'
        neumann_nodes(1, :) = (num_nodes - Nr):1:num_nodes;
        neumann_nodes(2,:) = right_value .* ones(1, Nr + 1);
    otherwise
        disp('You entered an incorrect field for the type of BC on the right boundary.');
end

a_k = [];

if isempty(dirichlet_nodes)
    disp('no dirichlet nodes')
else
    a_k = dirichlet_nodes(2,:)';
end

