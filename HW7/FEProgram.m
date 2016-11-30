clear all

L = 1.0;                        % problem domain
left = 'Dirichlet';             % left boundary condition 
left_value = 0.5;               % left Dirichlet boundary condition value
right = 'Neumann';              % right boundary condition type
right_value = -(5e-6);          % right Dirichlet boundary condition value
fontsize = 16;                  % fontsize for plots
num_elem = 100;                 % number of finite elements
shape_order = 2;                % linear elements
end_time = 6.5e3;               % end simulation time
num_steps = 1000;                % number of time steps
dt = end_time / num_steps;      % time step size
ic = 0.5;                       % initial condition 
discr = 10;                     % plot every discr time steps

% specify D and Tau over the domain in a block structure
D_blocks    = [2.4 2.0 1.5 0.6 1.3 0.14 1.1 2.2 2.0 1.5].* (10^(-6));
Tau_blocks  = [1.2 0.8 0.3 1.4 1.15 0.75 0.35 0.85 1.25 2.0].* (10^(-3));
space_blocks = 0.1:0.1:L;

% form the permutation matrix for assembling the global matrices
[permutation] = permutation(shape_order);    

parent_domain = -1:0.1:1;
physical_domain = linspace(0, L, num_elem * length(parent_domain) - (num_elem - 1));    

% define the quadrature rule
[wt, qp] = quadrature(shape_order);

% interpolate D and Tau into the physical domain
[D_physical_domain] = PhysicalInterpolation(physical_domain, space_blocks, D_blocks);
[Tau_physical_domain] = PhysicalInterpolation(physical_domain, space_blocks, Tau_blocks);

% perform the meshing
[num_nodes, num_nodes_per_element, LM, coordinates] = mesh(L, num_elem, shape_order);

% specify the boundary conditions
[dirichlet_nodes, neumann_nodes, a_k] = BCnodes(left, right, left_value, right_value, num_nodes);
num_dirichlet = length(dirichlet_nodes(1,:));

n = 1; % index for the time step

% apply the initial condition
soln_condensed_cell = cell([1, num_steps]);
soln_condensed_cell{1, n} = ic .* ones(1, num_nodes - num_dirichlet)';

% initialize the solution and derivative cells
soln_FE_cell = cell([1, length(physical_domain)]);
soln_derivative_FE_cell = cell([1, length(physical_domain)]);

soln_FE_cell{1, n} = ic .* ones(1, length(physical_domain));
soln_derivative_FE_cell{1, n} = zeros(1, length(physical_domain));

% interpolate D and Tau into the an elemental basis
[D_elem, right_endpoint_index, right_endpoint_coordinate] = ElementInterpolation(coordinates, num_elem, num_nodes_per_element, space_blocks, D_blocks);
[Tau_elem, right_endpoint_index, right_endpoint_coordinate] = ElementInterpolation(coordinates, num_elem, num_nodes_per_element, space_blocks, Tau_blocks);

K_cell = cell([1, num_elem]);
M_cell = cell([1, num_elem]);
F_cell = cell([1, num_elem]);

K = zeros(num_nodes, num_nodes);
M = zeros(num_nodes, num_nodes);
F = zeros(num_nodes, 1);

for elem = 1:num_elem
    k = zeros(num_nodes_per_element);
    m = zeros(num_nodes_per_element);
    f = zeros(num_nodes_per_element, 1);

     for l = 1:length(qp)
         for i = 1:num_nodes_per_element
             [N, dN, x_xe, dx_dxe] = shapefunctions(qp(l), shape_order, coordinates, LM, elem);

             % assemble the (elemental) forcing vector
             if strcmp(right, 'Neumann')
                 if ((neumann_nodes(1,1) == (elem + 1)) && (i == num_nodes_per_element) && (l == 1))
                    f(i) = f(i) - neumann_nodes(2, 1);
                 end
             end

             for j = 1:num_nodes_per_element
                 % assemble the (elemental) stiffness matrix
                 k(i,j) = k(i,j) + wt(l) * (D_elem(elem) * dN(i) * dN(j) / dx_dxe + Tau_elem(elem) * N(i) * N(j) * dx_dxe);
                 
                 % assemble the (elemental) mass matrix
                 m(i,j) = m(i,j) + wt(l) * N(i) * N(j) * dx_dxe;
             end
         end
     end

     % store elemental values into cells
     K_cell{1, elem} = k;
     M_cell{1, elem} = m;
     F_cell{1, elem} = f;   
end
    
% assemble into the global matrices
for elem = 1:num_elem
     m = 1;
     for m = 1:length(permutation(:,1))
        i = permutation(m,1);
        j = permutation(m,2);
        K(LM(elem, i), LM(elem, j)) = K_cell{1, elem}(i, j) + K(LM(elem, i), LM(elem, j));
        M(LM(elem, i), LM(elem, j)) = M_cell{1, elem}(i, j) + M(LM(elem, i), LM(elem, j));
     end
     
     for i = 1:length(f)
        F(LM(elem, i)) = F((LM(elem, i))) + F_cell{1,elem}(i);
     end 
end

% perform static condensation to remove known Dirichlet nodes from solve
[K_uu, K_uk, F_u, F_k] = condensation(K, F, num_nodes, dirichlet_nodes);
[M_uu, M_uk, F_u, F_k] = condensation(M, F, num_nodes, dirichlet_nodes);

for n = 1:num_steps
    % perform the very first solve using Gaussian elimination (time-dependent case)
    A_mat = (1 / dt) * M_uu + K_uu;
    b_mat = F_u + (1 / dt) * (M_uu * soln_condensed_cell{1,n} + M_uk * dirichlet_nodes(2,:)') - ((1 / dt) * M_uk + K_uk) * dirichlet_nodes(2,:)';
    a_u_condensed = A_mat \ b_mat;
    soln_condensed_cell{1, n + 1} = a_u_condensed;

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
    soln_FE_cell{n+1} = solution_FE;
    soln_derivative_FE_cell{n+1} = solution_derivative_FE;
    
    n = n + 1;
end

[none] = PlotInTime(1, soln_FE_cell, physical_domain, 'Solution', discr, num_steps, 1);
%[none] = PlotInTime(1, soln_derivative_FE_cell, physical_domain, 'Derivative of the Solution', discr, num_steps, 2);
[none] = PlotInTime(-D_physical_domain, soln_derivative_FE_cell, physical_domain, 'Flux = -D dc/dx', discr, num_steps, 2);
saveas(gcf, '1000', 'jpeg')
