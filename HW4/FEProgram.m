clear all

% select which type of plot you want to make - at least one flag must equal 1
N_plot_flag = 0;                    % 1 - plot the solutions for various N

L = 1.0;                            % problem domain
shape_order = 2;                    % number of nodes per element
E = 1.0;                            % elastic modulus
left = 'Dirichlet';                 % left BC 
left_value = 1.0;                   % left Dirichlet BC value
right = 'Dirichlet';                % right BC type
right_value = cos(10 * pi * L^5);   % right Dirichlet BC value
tolerance = 0.05;                   % convergence tolerance
refine_tol = 0.05;                  % refinement tolerance
energy_norm = tolerance + 1;        % arbitrary initialization value
fontsize = 16;                      % fontsize for plots
num_refinements = 0;                % number of refinements to prevent loop
max_refinements = 1;                % maximum number of refinements + 1

% form the permutation matrix for assembling the global matrices
[permutation] = permutation(shape_order);

N_elem = [20];

% index for collecting error
finished_refining = 0;
e = 1;

for num_elem = N_elem
% uncomment to find how many elements are required to reach the error
% tolerance
 %while energy_norm > tolerance
%     num_elem = num_elem + 1;
    
    % --- ANALYTICAL SOLUTION --- %
    parent_domain = -1:0.01:1;
    physical_domain = zeros(1, num_elem * length(parent_domain) - (num_elem - 1));
   
    % perform the meshing
    [num_nodes, num_nodes_per_element, LM, coordinates] = mesh(L, num_elem, shape_order);
    
    %coordinates(:,1) = [0, linspace(0.5, L, num_elem)];
    
    
    
    % --- ADAPTIVE MESH REFINEMENT --- %
    while ((finished_refining ~= 1) && (num_refinements <= max_refinements))
        
            % create a new physical domain
            j = 1;
            for elem = 1:num_elem
                discretization = linspace(coordinates(LM(elem, 1)), coordinates(LM(elem, num_nodes_per_element)), length(parent_domain));
                if elem == 1
                    physical_domain(j:length(parent_domain)) = discretization;
                    j = j + length(parent_domain);
                else
                    physical_domain(j:(j + length(parent_domain) - 2)) = discretization(2:end);
                    j = j + length(parent_domain) - 1;
                end

            end
    
            % --- ANALYTICAL SOLUTION --- %
            solution_analytical = cos(10 .* pi .* physical_domain .^ 5);
            solution_analytical_derivative = - 10 .* pi .* 5 .* physical_domain .^ 4 .* sin(10 .* pi .* physical_domain .^ 5);
        
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
                         f(i) = f(i) - wt(l) * - E * (200 * pi * x_xe ^3 * sin(10 * pi * x_xe ^ 5) + 2500 * pi^2 * x_xe^8 * cos(10 * pi * x_xe^5)) * N(i) * dx_dxe;

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
            
            % check the derivatives
            %plot(physical_domain, solution_derivative_FE, 'r*', physical_domain, solution_analytical_derivative, 'k')
            
            % compute the energy norm
            energy_norm_bottom = sqrt(trapz(physical_domain, solution_analytical_derivative .* E .* solution_analytical_derivative));
            energy_norm_top = sqrt(trapz(physical_domain, (solution_derivative_FE - solution_analytical_derivative) .* E .* (solution_derivative_FE - solution_analytical_derivative)));
            energy_norm = energy_norm_top ./ energy_norm_bottom;

            % determine the indices for the physical_domain to line up with coordinates
            k = 1;
            j = 1;
            bounds = zeros(1, length(coordinates));
            for i = 1:length(physical_domain)
                if (abs(coordinates(j) - physical_domain(i)) < 1e-6)
                    bounds(k) = i;
                    j = j + 1;
                    k = k + 1;
                end
            end

            if (N_plot_flag)
                plot(physical_domain, solution_FE)
                hold on
            end

            % uncomment to find out how many elements are needed to reach the error
            % tolerance
            % end
            e_N(e) = energy_norm;
            e = e + 1;

            if (N_plot_flag)
                plot(physical_domain, solution_analytical, 'k')
                txt = cell(length(N_elem),1);
                for i = 1:length(N_elem)
                   txt{i}= sprintf('N = %i', N_elem(i));
                end
                txt{i+1} = 'analytical';
                h = legend(txt);
                xlabel('Problem domain', 'FontSize', fontsize)
                ylabel('Solution', 'FontSize', fontsize)
                %saveas(gcf, 'Nplot', 'jpeg')
                close all
            end
            
            % compute the energy norm over each element
            eN_per_elem = zeros(1, length(num_elem)); % works only for linear elements
            elem_length = zeros(1, length(num_elem));
            for i = 1:num_elem
                elem_length(i) = coordinates(i+1, 1) - coordinates(i,1);
                spatial_domain = physical_domain(bounds(i):bounds(i+1));
                sprintf('Start: %.6f, End: %.6f', spatial_domain(1), spatial_domain(end));
                
                dFE = solution_derivative_FE(bounds(i):bounds(i+1));
                dAN = solution_analytical_derivative(bounds(i):bounds(i+1));
                
                eN_per_elem(i) = trapz(spatial_domain, (dFE - dAN) .* E .* (dFE - dAN));
            end
            
            sprintf('Difference between sum and exact = %.6f', sum(eN_per_elem) - energy_norm_top^2)
            
            % plot A_I as a function of the element number
            A_I = sqrt((1 ./ elem_length) .* eN_per_elem ./ ((1 ./ L) .* energy_norm_bottom .^ 2));
            coords = coordinates(:,1);
            plot(coords(2:1:end), A_I, '*-', physical_domain, solution_analytical, 'k')
            hold on
            xlabel('Element Number' , 'FontSize', fontsize)
            ylabel(sprintf('A_I for %i Elements', num_elem), 'FontSize', fontsize)
            %saveas(gcf, 'A_I_NoRefinement', 'jpeg')

            % determine which elements need to be refined
            clearvars refine
            j = 1;
            for i = 1:length(A_I)
                if (A_I(i) < refine_tol)
                    % no refinement
                else
                    % refinement
                    refine(j) = i;
                    j = j + 1;
                end
            end
                   
            if (j == 1)
                finished_refining = 1;
            else 
                num_refinements = num_refinements + 1;
            
                % update coordinates vector (only works for linear elements)
                num_elem_new = num_elem + length(refine);
                num_nodes_new = (shape_order - 1) * num_elem_new + 1;
                coordinates_new = zeros(num_nodes_new, 3);

                j = 1;
                k = 1;
                l = 1;
                for i = 1:num_elem
                    if ((j <= length(refine)) && (refine(j) == i)) % refine this element
                        increment = 0.5 * (coordinates(k+1,1) - coordinates(k,1));
                        coordinates_new(l,1) = coordinates(k,1);
                        coordinates_new(l+1,1) = coordinates(k,1) + increment;
                        coordinates_new(l+2,1) = coordinates(k+1,1);
                        k = k + 1;
                        l = l + 2;
                        j = j + 1;
                    else
                        coordinates_new(l,1) = coordinates(k,1);
                        coordinates_new(l+1,1) = coordinates(k+1,1);
                        l = l + 1;
                        k = k + 1;
                    end
                end

                coordinates = coordinates_new;
                num_elem = num_elem_new;
                num_nodes = num_nodes_new;

                % update LM
                num_nodes_per_element = shape_order;

                LM = zeros(num_elem, num_nodes_per_element); 

                for i = 1:num_elem
                    for j = 1:num_nodes_per_element
                        LM(i,j) = num_nodes_per_element * (i - 1) + j - (i - 1);
                    end
                end
        end
    end
end

% uncomment to find out how many elements are needed to reach the error
% tolerance
%sprintf('Number elements needed: %i', num_elem)