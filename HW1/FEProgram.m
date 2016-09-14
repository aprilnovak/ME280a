clear all

L = 1.0;                % problem domain
k_freq = 1.0;           % forcing frequency
num_elem = 10.0;         % number of finite elements
shape_order = 3;        % number of nodes per element
E = 1.0;                % elastic modulus
left = 'Dirichlet';     % left boundary condition 
left_value = 0.0;       % left Dirichlet boundary condition value
right = 'Dirichlet';    % right boundary condition type
right_value = 1.0;      % right Dirichlet boundary condition value

% perform the meshing
[num_nodes, num_nodes_per_element, LM, coordinates] = mesh(L, num_elem, shape_order);

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

% apply the boundary conditions
[dirichlet_nodes, neumann_nodes] = BCnodes(left, right, num_nodes);

% perform static condensation to remove known Dirichlet nodes from solve
[K_condensed, F_condensed] = condensation(K, F, num_nodes, dirichlet_nodes);

% perform the solve
a_condensed = K_condensed \ F_condensed;

% expand a_condensed to include the Dirichlet nodes






