clear all

L = 1.0;                        % problem domain
k_freq = 4.0;                   % forcing frequency
num_elem = 10;                  % number of finite elements (initial guess)
shape_order = 3;                % number of nodes per element
E = 0.1;                        % elastic modulus
left = 'Dirichlet';             % left boundary condition 
left_value = 0.0;               % left Dirichlet boundary condition value
right = 'Dirichlet';            % right boundary condition type
right_value = 1.0;              % right Dirichlet boundary condition value
tolerance = 0.05;               % convergence tolerance
energy_norm = tolerance + 1;    % arbitrary initialization value

% form the permutation matrix for assembling the global matrices
[permutation] = permutation(shape_order);

% while energy_norm > tolerance
%     if (energy_norm - tolerance < 0.005)
%         num_elem = num_elem + 1     % increase by 1 element if very near the convgergence
%     elseif (energy_norm - tolerance < 0.01)
%         num_elem = num_elem + 5     % increase by 5 elements if near the convergence
%     else
%         num_elem = num_elem + 10    % increase by 10 elements each time for a faster run
%     end

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

%     syms xe
% 
%     switch shape_order
%         case 2
%             N1(xe) = (1 - xe) ./ 2;
%             N2(xe) = (1 + xe) ./ 2;
%             dN1 = - 1/2;
%             dN2 = 1/2;
%             N = {N1, N2};
%             dN = {dN1, dN2};
%         case 3
%             N1(xe) = xe .* (xe - 1) ./ 2;
%             N2(xe) = - (xe - 1) .* (1 + xe);
%             N3(xe) = xe .* (1 + xe) ./ 2;
%             dN1(xe) = xe - 1/2;
%             dN2(xe) = -2 .* xe;
%             dN3(xe) = 1/2 + xe;
%             N = {N1, N2, N3};
%             dN = {dN1, dN2, dN3};
%         otherwise
%             disp('You entered an unsupported shape function order.');
%     end
% 
%     i = 1;
%     for elem = 1:num_elem
% 
%         if (shape_order == 2)
%             r(xe) = coordinates(LM(elem, 1), 1)*N{1}(xe) + coordinates(LM(elem, 2), 1)*N{2}(xe);
%             J = coordinates(LM(elem, 1), 1) * dN{1} + coordinates(LM(elem, 2), 1) * dN{2};
%             solution = a(LM(elem, 1))*N{1}(parent_domain) + a(LM(elem, 2))*N{2}(parent_domain);
%             solution_derivative = a(LM(elem, 1))*dN{1} + a(LM(elem, 2))*dN{2};
%         else
%             r(xe) = coordinates(LM(elem, 1), 1)*N{1}(xe) + coordinates(LM(elem, 2), 1)*N{2}(xe) + coordinates(LM(elem, 3), 1)*N{3}(xe);
%             J = coordinates(LM(elem, 1), 1) * dN{1}(xe) + coordinates(LM(elem, 2), 1) * dN{2}(xe) + coordinates(LM(elem, 3), 1) * dN{3}(xe);
%             solution = a(LM(elem, 1))*N{1}(parent_domain) + a(LM(elem, 2))*N{2}(parent_domain) + a(LM(elem, 3))*N{3}(parent_domain);
%             solution_derivative = a(LM(elem, 1))*dN{1}(parent_domain) + a(LM(elem, 2))*dN{2}(parent_domain) + a(LM(elem, 3))*dN{3}(parent_domain);
%         end
% 
%         % sample the solution into a vector u_sampled_solution
%         u_sampled_solution_matrix(i,:) = solution;
%         u_sampled_solution_derivative_matrix(i,:) = solution_derivative ./ J;
%         i = i + 1;
% 
% %         if elem == 1
% %             plot(coordinates(:,1), a, 'ro','LineWidth',2)
% %             grid on
% %             hold on
% %         end
% % 
% %         plot(r(parent_domain), solution, '-k')
% %         hold on
%     end
% 
%     %plot the analytical solution
% %     plot(physical_domain, solution_analytical)
% %     hold on
% 
%     % assemble u_sampled_solution into a single vector
%     j = length(u_sampled_solution_matrix(1,:)) + 1;
% 
%     for i = 1:length(u_sampled_solution_matrix(:,1))
%         if i == 1
%             solution_FE(1:length(u_sampled_solution_matrix(i,:))) = u_sampled_solution_matrix(i,:);
%             if shape_order == 2
%                 solution_derivative_FE(1:length(u_sampled_solution_matrix(i,:))) = u_sampled_solution_derivative_matrix(i,:) * ones(1,length(u_sampled_solution_matrix(i,:)));
%             else
%                 solution_derivative_FE(1:length(u_sampled_solution_derivative_matrix(i,:))) = u_sampled_solution_derivative_matrix(i,:);
%             end
%         else
%             solution_FE(j:(j + length(u_sampled_solution_matrix(1,:)) - 2)) = u_sampled_solution_matrix(i,2:end);
%             if shape_order == 2
%                 solution_derivative_FE(j:(j + length(u_sampled_solution_matrix(1,:)) - 2)) = u_sampled_solution_derivative_matrix(i,:) * ones(1,length(u_sampled_solution_matrix(i,2:end)));
%             else
%                 solution_derivative_FE(j:(j + length(u_sampled_solution_derivative_matrix(1,:)) - 2)) = u_sampled_solution_derivative_matrix(i,2:end);
%             end
%             j = j + length(u_sampled_solution_matrix(1,:)) - 1;
%         end
%     end
% % 
% %     plot(physical_domain, solution_derivative_FE, 'go')
% %     hold on
% %     plot(physical_domain, solution_analytical_derivative, 'm-')
% 
%     % compute the energy norm
%     energy_norm_analytical = trapz(physical_domain, solution_analytical_derivative .* E .* solution_analytical_derivative);
%     energy_norm_FE = trapz(physical_domain, solution_derivative_FE .* E .* solution_derivative_FE);
%     % change to eval for higher-order elements
%     energy_norm = sqrt(energy_norm_analytical - energy_norm_FE) ./ sqrt(energy_norm_analytical);
% %end

%[dummy] = FEplot(a, coordinates, shape_order, num_elem, physical_domain, LM, num_nodes_per_element);

b = zeros(1, shape_order);
A = zeros(shape_order);
m = length(parent_domain) + 1;
p = 1;

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
    u_sampled_solution_derivative_matrix2(p,:) = derivative_over_element;
    p = p + 1;
    
    
%     plot(element_domain, derivative_over_element, 'k*')
%     hold on
       
end

% assemble solution into a single vector (for plotting and computing the
% energy norm)

for i = 1:length(u_sampled_solution_matrix(:,1))
    if i == 1
        solution_FE(1:length(u_sampled_solution_matrix(i,:))) = u_sampled_solution_matrix(i,:);
        if shape_order == 2
            solution_derivative_FE(1:length(u_sampled_solution_matrix(i,:))) = u_sampled_solution_derivative_matrix2(i,:);
        else
            solution_derivative_FE(1:length(u_sampled_solution_derivative_matrix2(i,:))) = u_sampled_solution_derivative_matrix2(i,:);
        end
    else
        solution_FE(m:(m + length(u_sampled_solution_matrix(1,:)) - 2)) = u_sampled_solution_matrix(i,2:end);
        if shape_order == 2
            solution_derivative_FE(m:(m + length(u_sampled_solution_matrix(1,:)) - 2)) = u_sampled_solution_derivative_matrix2(i,2:end);
        else
            solution_derivative_FE(m:(m + length(u_sampled_solution_derivative_matrix2(1,:)) - 2)) = u_sampled_solution_derivative_matrix2(i,2:end);
        end
        m = m + length(u_sampled_solution_matrix(1,:)) - 1;
    end
end

plot(physical_domain, solution_FE, 'r')
hold on
plot(physical_domain, solution_analytical, 'k')

