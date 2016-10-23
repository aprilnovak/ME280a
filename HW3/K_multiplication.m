function [result] = K_multiplication(K, a, dirichlet_nodes)
% This function performs multiplication of a matrix K with a vector a

% result = zeros(1, length(K(:,1)));
% for i = 1:length(K(:,1))
%     for j = 1:length(K(1,:))
%     result(i) = result(i) + K(i,j) * a(j);
%     end
% end
% 
% result = result';



K_row = 1;
result = zeros(1, length(K(:,1)));

% do the multiplication
for i = 1:length(K(:,1)) % row index
    if (find(dirichlet_nodes(1,:) == i))
    else
        for j = 1:length(K(1,:)) % column index
            if (find(dirichlet_nodes(1,:) == j))
            else
                result(i) = result(i) + K(i,j) * a(j);
            end
        end
    end
end

% remove the dirichlet node positions from result
new_result = ones(length(K(1,:)) - length(dirichlet_nodes(1,:)), 1);

m = 1;
for i = 1:length(result)
    if (find(dirichlet_nodes(1,:) == i))
    else
        new_result(m) = result(i);
        m = m + 1;
    end
end

result = new_result;

end

