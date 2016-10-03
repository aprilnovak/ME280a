clear all

% select which type of plot you want to make - at least one flag must equal 1
k_plot_flag = 1;                % 1 - plot the error as a function of N
N_plot_flag = 0;                % 1 - plot the solutions for various N

L = 1.0;                        % problem domain
k_freq = 2;                     % forcing frequency
num_elem = 5;                   % number of finite elements (initial guess)
shape_order = 4;                % number of nodes per element
E = 0.2;                        % elastic modulus
left = 'Dirichlet';             % left boundary condition 
left_value = 3.0;               % left Dirichlet boundary condition value
right = 'Dirichlet';            % right boundary condition type
right_value = -1.0;              % right Dirichlet boundary condition value
tolerance = 0.05;               % convergence tolerance
energy_norm = tolerance + 1;    % arbitrary initialization value
fontsize = 16;                  % fontsize for plots

% form the permutation matrix for assembling the global matrices
[permutation] = permutation(shape_order);

if (N_plot_flag)
     N_elem = [100];
elseif (k_plot_flag)
     N_elem = 1:500;
else
    disp('Either N_plot_flag or k_plot_flag has to equal 1.');
end

for k_freq = 12

% index for collecting error
e = 1;

for num_elem = N_elem
% uncomment to find how many elements are required to reach the error
% tolerance
% while energy_norm > tolerance
%     num_elem = num_elem + 1;
    
    % --- ANALYTICAL SOLUTION --- %
    parent_domain = -1:0.01:1;
    physical_domain = linspace(0, L, num_elem * length(parent_domain) - (num_elem - 1));
    gamma = 2 * pi * k_freq ./ L;
    term1 = 2 * (k_freq .^ 3) * sin(gamma .* physical_domain) ./ (E .* gamma .^3);
    term2 = (k_freq .^ 3) * physical_domain .* cos(gamma .* physical_domain) ./ (E .* gamma .^ 2);
    C_1 = left_value;
    C_2 = (right_value - C_1 - (2 * (k_freq .^ 3) * sin(gamma .* L) ./ (E .* gamma .^3)) + ((k_freq .^ 3) * L .* cos(gamma .* L) ./ (E .* gamma .^ 2))) ./ L;
    solution_analytical = C_1 + term1 - term2 + C_2 .* physical_domain;
    term1_1 = 2 * (k_freq .^ 3) * cos(gamma .* physical_domain) .* gamma ./ (E .* gamma .^3);
    term2_1 = ((k_freq .^ 3) ./ (E .* gamma .^ 2)) .* (physical_domain .* gamma .* - sin(gamma .* physical_domain) + cos(gamma .* physical_domain));
    solution_analytical_derivative = (k_freq^3) * (gamma .* physical_domain .* sin(gamma .* physical_domain) + cos(gamma .* physical_domain)) ./ (E .* gamma .^ 2);
    
    
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
                 f(i) = f(i) - wt(l) * x_xe * (k_freq .^ 3) * cos(gamma * x_xe) * N(i) * dx_dxe;

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

% assemble the solution in the physical domain
[solution_FE, solution_derivative_FE] = postprocess(num_elem, parent_domain, a, LM, num_nodes_per_element, shape_order, coordinates, physical_domain);

% compute the energy norm
energy_norm_bottom = sqrt(trapz(physical_domain, solution_analytical_derivative .* E .* solution_analytical_derivative));
energy_norm_top = sqrt(trapz(physical_domain, (solution_derivative_FE - solution_analytical_derivative) .* E .* (solution_derivative_FE - solution_analytical_derivative)));
energy_norm = energy_norm_top ./ energy_norm_bottom;
sprintf('energy norm: %f', energy_norm)

if (N_plot_flag)
    plot(physical_domain, solution_FE)
    hold on
end

% uncomment to find out how many elements are needed to reach the error
% tolerance
%end
e_N(e) = energy_norm;
e = e + 1;

end

if (N_plot_flag)
    plot(physical_domain, solution_analytical)
    h = legend('N = 2', 'N = 4', 'N = 8', 'N = 16', 'N = 32', 'N = 64', 'N = 128','analytical', 'Location', 'southeast');
    set(h, 'FontSize', fontsize - 2);
    xlabel('Problem domain', 'FontSize', fontsize)
    ylabel(sprintf('Solution for k = %i', k_freq), 'FontSize', fontsize)
    text(0.85, 0.5, sprintf('k = %i', k_freq), 'FontSize', fontsize, 'FontWeight', 'bold', 'EdgeColor', [0 0 0])
    saveas(gcf, sprintf('Nplot_for_k_%i', k_freq), 'jpeg')
    %close all
end

if (k_plot_flag)
    loglog(N_elem, e_N, '*-')
    hold on
    xlabel('Number of elements', 'FontSize', fontsize)
    ylabel('Energy norm', 'FontSize', fontsize)
end

end

if (k_plot_flag)
    h2 = legend('k = 1', 'k = 2', 'k = 4', 'k = 8', 'k = 16', 'k = 32');
    set(h2, 'FontSize', fontsize);
    saveas(gcf, 'eN_vs_N', 'jpeg')
end


% uncomment to find out how many elements are needed to reach the error
% tolerance
%sprintf('For k = %i, number elements: %i', k_freq, num_elem)