function [result] = K_multiplication(K, vector, dirichlet_nodes, non_dn)
% This function performs multiplication of a matrix K with a vector
result = zeros(1, length(K(:,1)));

for i = 1:length(K(:,1)) % row index
    for j = non_dn % column index
        result(i) = result(i) + K(i,j) * vector(j);
    end
end

result = CutoffDirichlet(result, dirichlet_nodes);
end

