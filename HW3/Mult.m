function [result] = Mult(K_cell, num_elem, vector, num_nodes, num_nodes_per_element, dirichlet_nodes, LM)
% This function multiplies the global K by a vector without needing to
% explicitly form the global K and by avoiding unnecessary multiplication
% of zeros. You don't need to cut out the Dirichlet nodes.

result = zeros(1, num_nodes);
for elem = 1:num_elem
    for i = 1:num_nodes_per_element
    	for j = 1:num_nodes_per_element
            %if find(dirichlet_nodes(1,:) == LM(elem, j))
            if ((LM(elem, j) == 1) || (LM(elem, j) == num_nodes)) % works good enough for only 2 dirichlet points
            else
                result(LM(elem,i)) = K_cell{1,elem}(i,j) * vector(LM(elem, j)) + result(LM(elem,i));
            end
        end
    end
end

result = CutoffDirichlet(result, dirichlet_nodes);

end

