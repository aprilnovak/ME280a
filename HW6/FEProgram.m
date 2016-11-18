clear all

L = pi;                         % problem domain (theta)
k_th = 2;                       % thermal conductivity
shape_order = 2;                % number of nodes per element
E = 0.1;                        % elastic modulus
To = 110;                       % temperature at theta = pi
left = 'Dirichlet';             % boundary condition at theta = pi
left_value = To;                % 
right = 'Neumann';              % boundary condition at theta = 0
right_value = 1.0;              % 
fontsize = 16;                  % fontsize for plots
Nr = 3;                         % number of radial layers
No = 12;                        % number of theta layers
N_elem = Nr * No;               % number of elements
num_nodes = (Nr + 1) * (No + 1);% number of nodes
ri = 3;                         % inner radius of arch
ro = 4;                         % outer radius of arch
dt = (ro - ri)/Nr;              % thickness of each radial layer

% form the permutation matrix for assembling the global matrices
[permutation] = permutation(shape_order);

for num_elem = N_elem
    
    % --- ANALYTICAL SOLUTION --- %
    parent_domain = -1:0.01:1;
    physical_domain = linspace(0, L, num_elem * length(parent_domain) - (num_elem - 1));
    C_o = 40 / k_th;
    C_1 = To - C_o * pi;
    solution_analytical = 10 .* sin(2 .* physical_domain) ./ k_th + C_o .* physical_domain + C_1;
    
    % for a 2-D mesh polar mesh
    [coordinates_polar, LM_polar] = polar_mesh(No, Nr, dt, num_nodes, ri, ro, num_elem);
    
%     % mesh of the semi-circle, with node numbering
%     x = coordinates_polar(:,1);
%     y = coordinates_polar(:,2);
%     scatter(x,y)
%     a = [1:num_nodes]';
%     b = num2str(a);
%     c = cellstr(b);
%     dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
%     text(x+dx, y+dy, c);
%     xlim([-ro-2*dx, ro+5*dx])
%     ylim([-ro/2, ro + ro/2])
%     %saveas(gcf, 'Mesh', 'jpeg')
    
%     % output the LM for report
%     for e = 1:length(LM_polar)
%        fprintf('%i & %i & %i & %i \\\\\n', LM_polar(e, 1), LM_polar(e,2), LM_polar(e,3), LM_polar(e,4))
%     end
    
    % original meshing (Cartesian)
    [num_nodes, num_nodes_per_element, LM, coordinates] = mesh(L, num_elem, shape_order);
    
    % specify the boundary conditions
    [dirichlet_nodes, neumann_nodes, a_k] = BCnodes(left, right, left_value, right_value, num_nodes, Nr, No);

    % define the quadrature rule
    [wt, qp] = quadrature(shape_order);

    % assemble the elemental k and elemental f
    K = zeros(num_nodes);
    F = zeros(num_nodes, 1);

    for elem = 1:num_elem
        k = zeros(num_nodes_per_element);
        f = zeros(num_nodes_per_element, 1);

        for ll = 1:length(qp) % eta loop
             for l = 1:length(qp) % xe loop
                 for i = 1:num_nodes_per_element
                     [N, dN_dxe, dN_deta, x_xe, dx_dxe] = shapefunctions(qp(l), qp(ll), shape_order, coordinates, LM, elem);

                     % assemble the (elemental) forcing vector
                     f(i) = f(i) - wt(l) * k_th * k_th * sin(2 * pi * k_th * x_xe / L) * N(i) * dx_dxe;

                     for j = 1:num_nodes_per_element
                         % assemble the (elemental) stiffness matrix
                         k(i,j) = k(i,j) + wt(l) * E * dN_dxe(i) * dN_dxe(j) / dx_dxe;
                     end
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

% plot(physical_domain, solution_FE)
% hold on



end

% plot(physical_domain, solution_analytical)
% h = legend('N = what','analytical', 'Location', 'southeast');
% set(h, 'FontSize', fontsize - 2);
% xlabel('Problem domain', 'FontSize', fontsize)
% ylabel(sprintf('Solution for k = %i', k_th), 'FontSize', fontsize)
%saveas(gcf, sprintf('Nplot_for_k_%i', k_th), 'jpeg')
%close all

