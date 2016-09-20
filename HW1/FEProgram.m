clear all

L = 2.0;                % problem domain
k_freq = 2.0;           % forcing frequency
num_elem = 4.0;         % number of finite elements
shape_order = 3;        % number of nodes per element
E = 0.1;                % elastic modulus
left = 'Dirichlet';     % left boundary condition 
left_value = 0.0;       % left Dirichlet boundary condition value
right = 'Dirichlet';    % right boundary condition type
right_value = 1.0;      % right Dirichlet boundary condition value

% perform the meshing
[num_nodes, num_nodes_per_element, LM, coordinates] = mesh(L, num_elem, shape_order);

% specify the boundary conditions
[dirichlet_nodes, neumann_nodes, a_k] = BCnodes(left, right, left_value, right_value, num_nodes);

% form the permutation matrix for assembling the global matrices
[permutation] = permutation(num_nodes_per_element);

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

% plot the solution in the parent domain, and then transform it to the
% physical domain

parent_domain = -1:0.1:1;

syms xe

switch shape_order
    case 2
        N1(xe) = (1 - xe) ./ 2;
        N2(xe) = (1 + xe) ./ 2;
        dN1(xe) = - 1/2;
        dN2(xe) = 1/2;
        N = {N1, N2};
        dN = {dN1, dN2};
    case 3
        N1(xe) = xe .* (xe - 1) ./ 2;
        N2(xe) = - (xe - 1) .* (1 + xe);
        N3(xe) = xe .* (1 + xe) ./ 2;
        dN1(xe) = xe - 1/2;
        dN2(xe) = -2 .* xe;
        dN3(xe) = 1/2 + xe;
        N = {N1, N2, N3};
        dN = {dN1, dN2, dN3};
    otherwise
        disp('You entered an unsupported shape function order.');
end

i = 1;
for elem = 1:num_elem
    r(xe) = coordinates(LM(elem, 1), 1)*N{1}(xe) + coordinates(LM(elem, 2), 1)*N{2}(xe) + coordinates(LM(elem, 3), 1)*N{3}(xe);
    J = diff(r,xe);
    solution = a(LM(elem, 1))*N{1}(parent_domain) + a(LM(elem, 2))*N{2}(parent_domain) + a(LM(elem, 3))*N{3}(parent_domain);
    
    % sample the solution into a vector u_sampled_solution
    u_sampled_solution_matrix(i,:) = solution;
    i = i + 1;
    
    if elem == 1
        plot(coordinates(:,1), a, 'ro','LineWidth',2)
        grid on
        hold on
    end
    
    plot(r(parent_domain), solution, '-k')
    hold on
end

% plot the analytical solution
physical_domain = linspace(0, L, num_elem * length(parent_domain) - (num_elem - 1));
C_1 = (right_value + (k_freq^2 * sin(2 * pi * k_freq) * (L / (2 * pi * k_freq))^2 / E)) / L;
solution_analytical = (1 ./ E) .* -k_freq.^2 .* sin(2 .* pi .* k_freq .* physical_domain ./ L) .* (L ./ (2 .* pi .* k_freq)).^2 + C_1 .* physical_domain + left_value;
plot(physical_domain, solution_analytical)

% assemble u_sampled_solution into a single vector
j = length(u_sampled_solution_matrix(1,:)) + 1;

for i = 1:length(u_sampled_solution_matrix(:,1))
    if i == 1 % first row
        solution_FE(1:length(u_sampled_solution_matrix(i,:))) = u_sampled_solution_matrix(i,:);
    else
        solution_FE(j:(j + length(u_sampled_solution_matrix(1,:)) - 2)) = u_sampled_solution_matrix(i,2:end);
        j = j + length(u_sampled_solution_matrix(1,:)) - 1;
    end
end

solution_difference = solution_FE - solution_analytical;


