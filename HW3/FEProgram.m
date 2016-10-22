clear all

% select which type of plot you want to make - at least one flag must equal 1
k_plot_flag = 0;                % 1 - plot the error as a function of order
k_plot_flag_dof = 0;            % 1 - plot the error as a function of DOF
N_plot_flag = 1;                % 1 - plot the solutions for various N

L = 1.0;                        % problem domain
k_freq = 12;                    % forcing frequency
num_elem = 5;                   % number of finite elements (initial guess)
shape_order = 2;                % number of nodes per element
E = 0.01;                       % elastic modulus
left = 'Dirichlet';             % left boundary condition 
left_value = -0.1;              % left Dirichlet boundary condition value
right = 'Dirichlet';            % right boundary condition type
right_value = 0.7;              % right Dirichlet boundary condition value
tolerance = 0.04;               % convergence tolerance
energy_norm = tolerance + 1;    % arbitrary initialization value
fontsize = 16;                  % fontsize for plots
pcg_error_tol = 0.000001;       % error tolerance for PCG

if (N_plot_flag)
     N_elem = [100];              % num_elem to cycle through for soln plots
elseif (k_plot_flag || k_plot_flag_dof)
     N_elem = 50:10:1000;        % num_elem to cycle through for e_N vs. N
else
    disp('Either N_plot_flag or k_plot_flag has to equal 1.');
end

Order = [4];              % shape function (orders - 1) to cycle thru

% specify E over the domain in a block structure
E_blocks = [2.5, 1.0, 1.75, 1.25, 2.75, 3.75, 2.25, 0.75, 2.0, 1.0];
space_blocks = 0.1:0.1:L;

if (space_blocks(end) ~= L)
    disp('Incorrect specifications for space_blocks.')
end

for shape_order = Order
    clearvars permutation

% form the permutation matrix for assembling the global matrices
[permutation] = permutation(shape_order);    

% index for collecting error
e = 1;

% initial guess for determining number of elements to reach error tol
%N_elem = 100;

for num_elem = N_elem

% uncomment to find how many elements are required to reach the error tol %
% while energy_norm > tolerance
%     num_elem = num_elem + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parent_domain = -1:0.1:1;
    physical_domain = linspace(0, L, num_elem * length(parent_domain) - (num_elem - 1));    

    % define the quadrature rule
    [wt, qp] = quadrature(shape_order);
    
    % interpolate E into the physical domain
    [E_physical_domain] = PhysicalInterpolation(physical_domain, space_blocks, E_blocks);
 
    % --- ANALYTICAL SOLUTION --- %
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
    
    % perform the meshing
    [num_nodes, num_nodes_per_element, LM, coordinates] = mesh(L, num_elem, shape_order);
    
    % interpolate E into the an elemental basis
    [E_elem, right_endpoint_index, right_endpoint_coordinate] = ElementInterpolation(coordinates, num_elem, num_nodes_per_element, space_blocks, E_blocks);

    % specify the boundary conditions
    [dirichlet_nodes, neumann_nodes, a_k] = BCnodes(left, right, left_value, right_value, num_nodes);

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
                     k(i,j) = k(i,j) + wt(l) * E_elem(elem) * dN(i) * dN(j) / dx_dxe;
                 end
             end
         end
         
         % place the elemental k matrix into the global K matrix
         m = 1;
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
[a_u_condensed, a_u_condensed_ge] = PCG(K_uu, F_u, K_uk, dirichlet_nodes, pcg_error_tol);
%plot(linspace(0, L, length(K_uu(1,:))), a_u_condensed_ge - a_u_condensed)

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
%sprintf('energy norm: %f', energy_norm)

if (N_plot_flag)
    plot(physical_domain, solution_FE)
    hold on
end

% uncomment to find how many elements are needed to reach the error tol   %
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e_N(e) = energy_norm;
e = e + 1;

end

if (N_plot_flag)
    plot(physical_domain, solution_analytical, 'k')
    txt = cell(length(N_elem),1);
    for i = 1:length(N_elem)
       txt{i}= sprintf('N = %i', N_elem(i));
    end
    txt{i+1} = 'analytical';
    h = legend(txt);
    set(h, 'FontSize', fontsize - 2);
    xlabel('Problem domain', 'FontSize', fontsize)
    ylabel(sprintf('Solution for order = %i', shape_order - 1), 'FontSize', fontsize)
    
    saveas(gcf, sprintf('Nplot_for_order_%i', shape_order - 1), 'jpeg')
    %close all
end

if (k_plot_flag || k_plot_flag_dof)
    if (k_plot_flag)
        independent_var = L ./ N_elem;
        independent_var_str = 'Element size h';
        filename = 'eN_vs_h';
    else
        independent_var = shape_order * N_elem;
        independent_var_str = 'Degrees of Freedom';
        filename = 'eN_vs_dof';
    end
    
    loglog(independent_var, e_N, '*-')
    hold on
    xlabel(independent_var_str, 'FontSize', fontsize)
    ylabel('Energy norm', 'FontSize', fontsize)
end

end

if (k_plot_flag || k_plot_flag_dof)
    txt = cell(length(Order),1);
    for i = 1:(length(Order))
        txt{i} = sprintf('Order = %i, |(slope)| = %i', Order(i) - 1, Order(i) - 1);
    end
    h2 = legend(txt);
    set(h2, 'FontSize', fontsize);
    saveas(gcf, filename, 'jpeg')
end

% uncomment to find out how many elements are needed to reach the error tol
%sprintf('For order = %i, number elements: %i', shape_order - 1, num_elem)