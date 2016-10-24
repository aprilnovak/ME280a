function [result] = CutoffDirichlet(vector, dirichlet_nodes)
% This function removes Dirichlet nodes from a vector

j = 1;
for i = 1:length(vector)
    if (find(dirichlet_nodes(1,:) == i))
    else
        result(j) = vector(i);
        j = j + 1;
    end
end

result = result';

end

