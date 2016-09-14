% Script to return the node numbers associated with different types of
% boundary conditions

function [dirichlet_nodes, neumann_nodes] = BCnodes(left, right, num_nodes)

dirichlet_nodes = int16.empty(1,0);
neumann_nodes = int16.empty(1,0);

i = 1;
switch left
    case 'Dirichlet'
        dirichlet_nodes(i) = 1;
    case 'Neumann'
        neumann_nodes(i) = 1;
    otherwise
        disp('You entered an incorrect field for the type of BC on the left boundary.');
end


i = i + 1;
switch right
    case 'Dirichlet'
        dirichlet_nodes(i) = num_nodes;
    case 'Neumann'
        neumann_nodes(i) = num_nodes;
    otherwise
        disp('You entered an incorrect field for the type of BC on the right boundary.');
end