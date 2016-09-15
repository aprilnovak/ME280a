clear all

L = 1.0;                % problem domain
k_freq = 1.0;           % forcing frequency
num_elem = 10.0;         % number of finite elements
shape_order = 3;        % number of nodes per element
E = 0.1;                % elastic modulus
left = 'Dirichlet';     % left boundary condition 
left_value = 0.0;       % left Dirichlet boundary condition value
right = 'Dirichlet';    % right boundary condition type
right_value = 1.0;      % right Dirichlet boundary condition value

% perform the meshing
[num_nodes, num_nodes_per_element, LM, coordinates] = mesh(L, num_elem, shape_order);

% specify the boundary conditions
[dirichlet_nodes, neumann_nodes] = BCnodes(left, right, left_value, right_value, num_nodes);

% form the permutation matrix for assembling the global matrices
[permutation] = permutation(num_nodes_per_element);

% define the quadrature rule
wt = [1.0, 1.0];                    % weights
qp = [-sqrt(1/3), sqrt(1/3)];       % points

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
             f(i) = f(i) + wt(l) * k_freq * k_freq * sin(2 * pi * k_freq * x_xe / L) * N(i) / dx_dxe;
             
             for j = 1:num_nodes_per_element
                 % assemble the (elemental) stiffness matrix
                 k(i,j) = k(i,j) + wt(l) * E * dN(i) * dN(j) * dx_dxe;
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
[K_condensed, F_condensed] = condensation(K, F, num_nodes, dirichlet_nodes);

% perform the solve
a_condensed = K_condensed \ F_condensed;

% expand a_condensed to include the Dirichlet nodes
a_expanded = zeros(num_nodes, 1);

a_row = 1;
i = 1;      % index for dirichlet_nodes
j = 1;      % index for expanded row
        
for a_row = 1:num_nodes
    if (find(dirichlet_nodes(1, :) == a_row))
        a_expanded(a_row) = dirichlet_nodes(2,i);
        i = i + 1;
    else
        a_expanded(a_row) = a_condensed(j);
        j = j + 1;
    end
end

% plot the solution in the parent domain, and then transform it to the
% physical domain
parent_domain = -1:0.1:1;

for elem = 1:num_elem
    solution_in_element = zeros(length(parent_domain), 1);
    for j = 1:length(parent_domain)
        [N, dN, x_xe, dx_dxe] = shapefunctions(parent_domain(j), shape_order, coordinates, LM, elem);
        for i = 1:num_nodes_per_element
            solution_in_element = solution_in_element + a_expanded(LM(elem, i)) * N(i);
        end
    end
end

% plot the analytical solution
x = 0:0.1:L;
C_1 = (right_value + (k_freq^2 * sin(2 * pi * k_freq) * (L / (2 * pi * k_freq))^2 / E)) / L;
solution_analytical = (1 / E) * -k_freq^2 * sin(2 * pi * k_freq * x / L) * (L / (2 * pi * k_freq))^2 + C_1 * x + left_value;

plot(x, solution_analytical, 'k')
%plot(coordinates(:,1), a_expanded, 'b')




