function [result] = AddDirichlet(vector, dirichlet_nodes)
% This function adds zero values for the locations of Dirichlet nodes

j = 1;
for i = 1:(length(vector) + length(dirichlet_nodes(1,:)))
    if (find(dirichlet_nodes(1,:) == i))
        result(i) = 0;
    else
        result(i) = vector(j);
        j = j + 1;
    end
end

result = result';

end

