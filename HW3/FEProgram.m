clear all

% select which type of plot you want to make - at least one flag must equal 1
k_plot_flag = 0;                % 1 - plot the error as a function of order
k_plot_flag_dof = 0;            % 1 - plot the error as a function of DOF
N_plot_flag = 1;                % 1 - plot the solutions for various N

L = 1.0;                        % problem domain
k_freq = 12;                    % forcing frequency
num_elem = 5;                   % number of finite elements (initial guess)
shape_order = 2;                % number of nodes per element
E = 0.01;                        % elastic modulus
left = 'Dirichlet';             % left boundary condition 
left_value = -0.3;              % left Dirichlet boundary condition value
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

Order = [3];              % shape function (orders - 1) to cycle thru


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
    
    % specify E over the domain in equally-spaced intervals
    E_blocks = [0.01, 0.01];
    space_blocks = [0.5, 1.0];

    % --- ANALYTICAL SOLUTION --- % (over the physical domain)
    gamma = 2 * pi * k_freq ./ L;
    term1 = 2 * (k_freq .^ 3) * sin(gamma .* physical_domain) ./ (E .* gamma .^3);
    term2 = (k_freq .^ 3) * physical_domain .* cos(gamma .* physical_domain) ./ (E .* gamma .^ 2);
    C_1 = left_value;
    C_2 = (right_value - C_1 - (2 * (k_freq .^ 3) * sin(gamma .* L) ./ (E .* gamma .^3)) + ((k_freq .^ 3) * L .* cos(gamma .* L) ./ (E .* gamma .^ 2))) ./ L;
    solution_analytical = C_1 + term1 - term2 + C_2 .* physical_domain;
    term1_1 = 2 * (k_freq .^ 3) * cos(gamma .* physical_domain) .* gamma ./ (E .* gamma .^3);
    term2_1 = ((k_freq .^ 3) ./ (E .* gamma .^ 2)) .* (physical_domain .* gamma .* - sin(gamma .* physical_domain) + cos(gamma .* physical_domain));
    solution_analytical_derivative = (k_freq^3) * (gamma .* physical_domain .* sin(gamma .* physical_domain) + cos(gamma .* physical_domain)) ./ (E .* gamma .^ 2) + C_2;
    
    % loop over the domain to find C_2 for each block
    C_2 = zeros(1, length(E_blocks));
    % assumes the domain begins at x=0 (also only works for this specific problem)
    % moves from left to right
    for i = 1:length(E_blocks) 
        if i == 1 % first block
            C_2(i) = left_value;
        else
            C_2(i) = (1/E_blocks(i-1)-1/E_blocks(i)) .* (-space_blocks(i) .* (k_freq .^ 3) .* cos(gamma .* space_blocks(i))./(gamma .^ 2) + 2 * k_freq .^ 3 .* sin(gamma .* space_blocks(i)) ./ (gamma .^ 3)) + C_2(i-1);
        end
    end
    
    % compute C_1 for the last block
    C_1 = (right_value - (-L .* (k_freq .^ 3) .* cos(gamma .* L)./(gamma .^ 2) + 2 * k_freq .^ 3 .* sin(gamma .* L) ./ (gamma .^ 3))/E_blocks(end) - C_2(end))/L;
    C_1_physical_domain = C_1 .* ones(1, length(physical_domain));
    
    % assembly E vector in physical_domain for the analytical solution
    j = 1;
    E_physical_domain = zeros(1,length(physical_domain));
    C_2_physical_domain = zeros(1,length(physical_domain));
    for i = 1:length(physical_domain)
        if (physical_domain(i) <= space_blocks(j))
            E_physical_domain(i) = E_blocks(j);
            C_2_physical_domain(i) = C_2(j);
        else
            j = j + 1;
            E_physical_domain(i) = E_blocks(j);
            C_2_physical_domain(i) = C_2(j);
        end
    end
    
    term1 = 2 * (k_freq .^ 3) * sin(gamma .* physical_domain) ./ (E_physical_domain .* gamma .^3);
    term2 = (k_freq .^ 3) * physical_domain .* cos(gamma .* physical_domain) ./ (E_physical_domain .* gamma .^ 2);
    solution_analytical_blocks = C_2_physical_domain + term1 - term2 + C_1_physical_domain .* physical_domain;
    
    term1_1 = 2 * (k_freq .^ 3) * cos(gamma .* physical_domain) .* gamma ./ (E_physical_domain .* gamma .^3);
    term2_1 = ((k_freq .^ 3) ./ (E_physical_domain .* gamma .^ 2)) .* (physical_domain .* gamma .* - sin(gamma .* physical_domain) + cos(gamma .* physical_domain));
    solution_analytical_derivative_blocks = (k_freq^3) * (gamma .* physical_domain .* sin(gamma .* physical_domain) + cos(gamma .* physical_domain)) ./ (E_physical_domain .* gamma .^ 2) + C_1_physical_domain;
    
    
    % perform the meshing
    [num_nodes, num_nodes_per_element, LM, coordinates] = mesh(L, num_elem, shape_order);

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
                     k(i,j) = k(i,j) + wt(l) * E * dN(i) * dN(j) / dx_dxe;
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
sprintf('energy norm: %f', energy_norm)

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
    close all
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