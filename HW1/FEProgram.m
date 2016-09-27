clear all

L = 1.0;                        % problem domain
k_freq = 1;                     % forcing frequency
num_elem = 3;                   % number of finite elements (initial guess)
shape_order = 2;                % number of nodes per element
E = 0.1;                        % elastic modulus
left = 'Dirichlet';             % left boundary condition 
left_value = 0.0;               % left Dirichlet boundary condition value
right = 'Dirichlet';            % right boundary condition type
right_value = 1.0;              % right Dirichlet boundary condition value
tolerance = 0.05;               % convergence tolerance
energy_norm = tolerance + 1;    % arbitrary initialization value

% form the permutation matrix for assembling the global matrices
[permutation] = permutation(shape_order);

% index for collecting error
e = 1;

% range of N for the simulation
N_elem = [2, 4, 8, 16, 32, 64, 128];

%for num_elem = N_elem

while energy_norm > tolerance
    num_elem = num_elem + 1;
    
    % --- ANALYTICAL SOLUTION --- %
    parent_domain = -1:0.01:1;
    physical_domain = linspace(0, L, num_elem * length(parent_domain) - (num_elem - 1));
    C_1 = (right_value + (k_freq^2 * sin(2 * pi * k_freq) * (L / (2 * pi * k_freq))^2 / E)) / L;
    solution_analytical = (1 ./ E) .* -k_freq.^2 .* sin(2 .* pi .* k_freq .* physical_domain ./ L) .* (L ./ (2 .* pi .* k_freq)).^2 + C_1 .* physical_domain + left_value;
    solution_analytical_derivative = -(1 ./ E) * k_freq * k_freq * cos(2 * pi * k_freq * physical_domain ./ L) * L ./ (2 * pi * k_freq) + C_1;

    % perform the meshing
    [num_nodes, num_nodes_per_element, LM, coordinates] = mesh(L, num_elem, shape_order);

    % specify the boundary conditions
    [dirichlet_nodes, neumann_nodes, a_k] = BCnodes(left, right, left_value, right_value, num_nodes);

    % define the quadrature rule
    [wt, qp] = quadrature(shape_order);

    % assemble the elemental k and elemental f
    K = zeros(num_nodes);
    F = zeros(num_nodes, 1);

    for elem = 1:num_elem
        k = zeros(num_nodes_per_element);
        f = zeros(num_nodes_per_element, 1);

         for l = 1:length(qp)
             for i = 1:num_nodes_per_element
                 [N, dN, x_xe, dx_dxe] = shapefunctions(qp(l), shape_order, coordinates, LM, elem);

                 % assemble the (elemental) forcing vector
                 f(i) = f(i) - wt(l) * k_freq * k_freq * sin(2 * pi * k_freq * x_xe / L) * N(i) * dx_dxe;

                 for j = 1:num_nodes_per_element
                     % assemble the (elemental) stiffness matrix
                     k(i,j) = k(i,j) + wt(l) * E * dN(i) * dN(j) / dx_dxe;
                 end
             end
         end

         % place the elemental k matrix into the global K matrix
         for m = 1:length(permutation(:,1))
            i = permutation(m,1);
            j = permutation(m,2);
            K(LM(elem, i), LM(elem, j)) = K(LM(elem, i), LM(elem, j)) + k(i,j);
         end

         % place the elemental f matrix into the global F matrix
         for i = 1:length(f)
            F(LM(elem, i)) = F((LM(elem, i))) + f(i);
         end
    end

% perform static condensation to remove known Dirichlet nodes from solve
[K_uu, K_uk, F_u, F_k] = condensation(K, F, num_nodes, dirichlet_nodes);

% perform the solve
a_u_condensed = K_uu \ (F_u - K_uk * dirichlet_nodes(2,:)');

% expand a_condensed to include the Dirichlet nodes
a = zeros(num_nodes, 1);

a_row = 1;
i = 1;      % index for dirichlet_nodes
j = 1;      % index for expanded row

for a_row = 1:num_nodes
    if (find(dirichlet_nodes(1, :) == a_row))
        a(a_row) = dirichlet_nodes(2,i);
        i = i + 1;
    else
        a(a_row) = a_u_condensed(j);
        j = j + 1;
    end
end

b = zeros(1, shape_order);
A = zeros(shape_order);
m = length(parent_domain) + 1;
p = 1;


% --- PLOTTING --- %
for elem = 1:num_elem
    % over each element, figure out the polynomial by solving a linear
    % system, Ax = b, where A depends on the order of the shape functions
    for i = 1:num_nodes_per_element
        b(i) = a(LM(elem, i));
    end
    
    for j = 1:shape_order % loop over the rows of A
        coordinate = coordinates(LM(elem, j));
        for l = 1:shape_order % loop over the columns of A
            A(j,l) = coordinate .^ (l - 1);
        end  
    end
    
    % solve for the coefficients on the actual polynomial
    coefficients = A\(b');
    
    % determine the solution over the element
    solution_over_element = zeros(1, length(parent_domain));
    element_domain = linspace(coordinates(LM(elem, 1)), coordinates(LM(elem, num_nodes_per_element)), length(parent_domain));
    for i = 1:num_nodes_per_element
        solution_over_element = solution_over_element + coefficients(i) .* (element_domain .^ (i - 1));
    end
    
    % determine the derivative over the element
    derivative_over_element = zeros(1, length(parent_domain));
    for i = 2:num_nodes_per_element % the derivative of the constant is zero
        derivative_over_element = derivative_over_element + coefficients(i) .* (i - 1) .* (element_domain .^ (i - 2));
    end
    
    % put into a matrix
    u_sampled_solution_matrix(p,:) = solution_over_element;
    u_sampled_solution_derivative_matrix(p,:) = derivative_over_element;
    p = p + 1;    
end

% assemble solution and derivative into a single vector
for i = 1:length(u_sampled_solution_matrix(:,1))
    if i == 1
        solution_FE(1:length(u_sampled_solution_matrix(i,:))) = u_sampled_solution_matrix(i,:);
        solution_derivative_FE(1:length(u_sampled_solution_derivative_matrix(i,:))) = u_sampled_solution_derivative_matrix(i,:);
    else
        solution_FE(m:(m + length(u_sampled_solution_matrix(1,:)) - 2)) = u_sampled_solution_matrix(i,2:end);
        solution_derivative_FE(m:(m + length(u_sampled_solution_derivative_matrix(1,:)) - 2)) = u_sampled_solution_derivative_matrix(i,2:end);
        m = m + length(u_sampled_solution_matrix(1,:)) - 1;
    end
end

% compute the energy norm
energy_norm_bottom = sqrt(trapz(physical_domain, solution_analytical_derivative .* E .* solution_analytical_derivative));
energy_norm_top = sqrt(trapz(physical_domain, (solution_derivative_FE - solution_analytical_derivative) .* E .* (solution_derivative_FE - solution_analytical_derivative)));
energy_norm = energy_norm_top ./ energy_norm_bottom;
sprintf('energy norm: %f', energy_norm)
e_N(e) = energy_norm;
e = e + 1;

% plot(physical_domain, solution_FE)
% hold on

end
% plot(physical_domain, solution_analytical, 'k')
% legend('N = 2', 'N = 4', 'N = 8', 'N = 16', 'N = 32', 'N = 64', 'N = 128', 'analytical solution', 'Location', 'southeast')
% xlabel('Problem domain')
% ylabel(sprintf('solution for k = %i', k_freq))
% 
% 
% figure()
% plot(N_elem, e_N)
% legend(sprintf('k = %i', k_freq))

sprintf('number elements: %i', num_elem)
% plot(physical_domain, solution_FE, 'r')
% hold on
% plot(physical_domain, solution_analytical, 'k')