function [result] = K_multiplication(K, vector, dirichlet_nodes)
% This function performs multiplication of a matrix K with a vector
result = zeros(1, length(K(:,1)));

for i = 1:length(K(:,1)) % row index
    for j = 1:length(K(1,:)) % column index
        if (find(dirichlet_nodes(1,:) == j))
        else
            result(i) = result(i) + K(i,j) * vector(j);
        end
    end
end

result = CutoffDirichlet(result, dirichlet_nodes);
end

