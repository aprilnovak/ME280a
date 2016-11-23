function [solution_analytical, solution_analytical_derivative, gamma] = AnalyticalSolution(L, space_blocks, left_value, E_blocks, E_physical_domain, right_value, k_freq, physical_domain)
gamma = 2 * pi * k_freq ./ L;

C_unknowns = zeros(1, length(space_blocks) + 1);
coeff_matrix = zeros(length(space_blocks) + 1, length(space_blocks) + 1);
coeff_vector = zeros(length(space_blocks) + 1, 1);

% left boundary condition
coeff_matrix(1,1) = 1; 
coeff_vector(1) = left_value;

% right boundary condition
coeff_matrix(end, end) = L / E_blocks(end);
coeff_matrix(end, end - 1) = 1;
coeff_vector(end) = right_value - (1 / E_blocks(end)) * RHSFunction(gamma, k_freq, space_blocks, E_blocks, length(space_blocks));

% intermediate conditions
for i = 1:(length(space_blocks) - 1)
    coeff_matrix(i+1, i) = 1;
    coeff_matrix(i+1, i+1) = -1;
    coeff_matrix(i+1, end) = (space_blocks(i) / E_blocks(i) - space_blocks(i) / E_blocks(i+1));
    coeff_vector(i+1) = - RHSFunction(gamma, k_freq, space_blocks, E_blocks, i) * (1/E_blocks(i) -1/E_blocks(i+1));
end

coeff_soln = coeff_matrix\coeff_vector;

% assemble block-oriented E and C_2 into physical_domain structure
[C_2_physical_domain2] = PhysicalInterpolation(physical_domain, space_blocks, coeff_soln(1:end-1));
C_1_physical_domain2 = coeff_soln(end) .* ones(1, length(physical_domain));

term1_2 = 2 * (k_freq .^ 3) * sin(gamma .* physical_domain) ./ (E_physical_domain .* gamma .^3);
term2_2 = (k_freq .^ 3) * physical_domain .* cos(gamma .* physical_domain) ./ (E_physical_domain .* gamma .^ 2);
solution_analytical = C_2_physical_domain2 + term1_2 - term2_2 + C_1_physical_domain2 .* physical_domain ./ E_physical_domain;

term1_1 = 2 * (k_freq .^ 3) * cos(gamma .* physical_domain) .* gamma ./ (E_physical_domain .* gamma .^3);
term2_1 = ((k_freq .^ 3) ./ (E_physical_domain .* gamma .^ 2)) .* (physical_domain .* gamma .* - sin(gamma .* physical_domain) + cos(gamma .* physical_domain));
solution_analytical_derivative = (k_freq^3) * (gamma .* physical_domain .* sin(gamma .* physical_domain) + cos(gamma .* physical_domain)) ./ (E_physical_domain .* gamma .^ 2) + C_1_physical_domain2 ./ E_physical_domain;



function [RHS_value] = RHSFunction(gamma, k_freq, space_blocks, E_blocks, node_number)
% function to repeatedly evaluate \circled{A} in the document

loc = space_blocks(node_number);
RHS_value = - loc * (k_freq^3) * cos(gamma * loc) / (gamma ^2) + 2 * k_freq ^ 3 * sin(gamma * loc) / (gamma ^ 3);
end


end

