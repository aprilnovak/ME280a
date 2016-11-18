function [coefficients_mat] = postprocess(num_elem, parent_domain, a, LM, num_nodes_per_elem, shape_order, coordinates, physical_domain)

b = zeros(1, num_nodes_per_elem);
A = zeros(num_nodes_per_elem);
coefficients_mat = zeros(num_elem, 4);
%m = length(parent_domain) + 1;
p = 1;

%u_sampled_solution_matrix = zeros(num_elem, length(parent_domain));
%u_sampled_solution_derivative_matrix = zeros(num_elem, length(parent_domain));

for elem = 1:num_elem
    % over each element, figure out the polynomial by solving a linear
    % system, Ax = b, where A depends on the order of the shape functions
    for i = 1:num_nodes_per_elem
        b(i) = a(LM(elem, i));
    end
    
    for j = 1:num_nodes_per_elem % loop over the rows of A
        x_coordinate = coordinates(LM(elem, j), 1);
        y_coordinate = coordinates(LM(elem, j), 2);
        A(j,:) = [1, x_coordinate, y_coordinate, x_coordinate * y_coordinate];
    end
    
    % solve for the coefficients on the actual polynomial
    coefficients = A\(b');
    
    % save those coefficients into a matrix (each row is an element)
    coefficients_mat(p,:) = coefficients;
    p = p + 1;
    
    % determine the solution over the element
%     solution_over_element = zeros(1, length(parent_domain));
%     element_domain = linspace(coordinates(LM(elem, 1)), coordinates(LM(elem, num_nodes_per_elem)), length(parent_domain));
%     for i = 1:num_nodes_per_elem
%         solution_over_element = solution_over_element + coefficients(i) .* (element_domain .^ (i - 1));
%     end
    
    % determine the derivative over the element
%     derivative_over_element = zeros(1, length(parent_domain));
%     for i = 2:num_nodes_per_elem % the derivative of the constant is zero
%         derivative_over_element = derivative_over_element + coefficients(i) .* (i - 1) .* (element_domain .^ (i - 2));
%     end
    
    % put into a matrix
%     u_sampled_solution_matrix(p,:) = solution_over_element;
%     u_sampled_solution_derivative_matrix(p,:) = derivative_over_element;
%     p = p + 1;    
end

% assemble solution and derivative into a single vector
% solution_FE = zeros(1, length(physical_domain));
% solution_derivative_FE = zeros(1, length(physical_domain));
% for i = 1:length(u_sampled_solution_matrix(:,1))
%     if i == 1
%         solution_FE(1:length(u_sampled_solution_matrix(i,:))) = u_sampled_solution_matrix(i,:);
%         solution_derivative_FE(1:length(u_sampled_solution_derivative_matrix(i,:))) = u_sampled_solution_derivative_matrix(i,:);
%     else
%         solution_FE(m:(m + length(u_sampled_solution_matrix(1,:)) - 2)) = u_sampled_solution_matrix(i,2:end);
%         solution_derivative_FE(m:(m + length(u_sampled_solution_derivative_matrix(1,:)) - 2)) = u_sampled_solution_derivative_matrix(i,2:end);
%         m = m + length(u_sampled_solution_matrix(1,:)) - 1;
%     end
% end