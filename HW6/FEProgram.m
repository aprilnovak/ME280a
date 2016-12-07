clear all

L = pi;                         % problem domain (theta)
k_th = 1;                       % thermal conductivity
shape_order = 2;                % number of nodes per element
E = 0.1;                        % elastic modulus
To = 110;                       % temperature at theta = pi
left = 'Dirichlet';             % boundary condition at theta = pi
left_value = To;                % 
right = 'Neumann';              % boundary condition at theta = 0
right_value = 1.0;              % dummy
fontsize = 16;                  % fontsize for plots
Nr = 1;                         % number of radial layers
No = 16;                        % number of theta layers
N_elem = Nr * No;               % number of elements
num_nodes = (Nr + 1) * (No + 1);% number of nodes
num_nodes_per_elem = 4;         % linear elements
ri = 3;                         % inner radius of arch
ro = 4;                         % outer radius of arch
dt = (ro - ri)/Nr;              % thickness of each radial layer

% form the permutation matrix for assembling the global matrices
[permutation] = permutation(num_nodes_per_elem);

for num_elem = N_elem
    
    % --- ANALYTICAL SOLUTION --- %
    parent_domain = -1:0.01:1;
    physical_domain = linspace(0, L, num_elem * length(parent_domain) - (num_elem - 1));
    C_o = 0;
    C_1 = To - C_o * pi;
    
    % for a 2-D mesh polar mesh
    [coordinates, LM] = polar_mesh(No, Nr, dt, num_nodes, ri, ro, num_elem);
    %[plot] = mesh_plots(coordinates, num_nodes, ro, LM);
    
    % specify the boundary conditions
    [dirichlet_nodes, neumann_nodes, a_k] = BCnodes(left, right, left_value, right_value, num_nodes, Nr, No);

    % define the quadrature rule
    [wt, qp] = quadrature(shape_order);

    % assemble the elemental k and elemental f
    K = zeros(num_nodes);
    F = zeros(num_nodes, 1);

    for elem = 1:num_elem
        k = 0;
        f = 0;

        for ll = 1:length(qp) % eta loop
             for l = 1:length(qp) % xe loop
                     [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(qp(l), qp(ll), num_nodes_per_elem, coordinates, LM, elem);
                     F_mat = transpose([dx_dxe, dx_deta; dy_dxe, dy_deta]);
                     J = det(F_mat);
                     r = sqrt(x_xe_eta^2 + y_xe_eta^2);
                     theta = acos(x_xe_eta / r);
                     
                     % assemble the (elemental) forcing vector
                     f = f + wt(ll) * wt(l) * (40 * sin(2 * theta) / (r^2)) * transpose(N) * J;
                     
                     % assemble the (elemental) stiffness matrix - correct
                     k = k + wt(ll) * wt(l) * transpose(inv(F_mat) * B) * k_th * inv(F_mat) * B * J;      
             end
             
             % apply flux boundary conditions (xe is constant, so this is
             % outside of the xe loop). Only the last Nr elements have BCs.
             if (num_elem - Nr) <= elem
                 q_flux = 20 / r;

                 % in the physical domain
                 N_hat = [0, -1];
                 % in the master domain
                 n_hat = [1, 0]';
                 
                 % using both quadrature points (incorrect?) gives good
                 % results...
                 f = f + wt(ll) * q_flux * transpose(N) * J * (N_hat * transpose(inv(F_mat)) * n_hat);
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
[mat] = postprocess(num_elem, parent_domain, a, LM, num_nodes_per_elem, shape_order, coordinates, physical_domain);

% subplot(1,2,1)
% start_theta = pi - pi/No;
% end_theta = pi;
% for i = 1:length(mat(:,1))
%     [r, theta] = meshgrid(ri:0.1:ro, linspace(start_theta, end_theta, length(ri:0.1:ro)));
%     x = r .* cos(theta);
%     y = r .* sin(theta);
%     z = mat(i,1) + mat(i,2).*x + mat(i,3).*y + mat(i,4).*x.*y;
%     surf(x,y,z, 'EdgeColor', 'none')
%     hold on
%     
%     if (mod(i, Nr) == 0)
%         end_theta = start_theta;
%         start_theta = start_theta - pi/No;
%     end
% end

% analytical solution
subplot(1,2,1)
ylabel('FE Solution')
start_theta = pi - pi/No;
end_theta = pi;
for i = 1:length(mat(:,1))
    [r, theta] = meshgrid(ri:0.1:ro, linspace(start_theta, end_theta, length(ri:0.1:ro)));
    x = r .* cos(theta);
    y = r .* sin(theta);
    z = 10 * sin(2*theta)/k_th + C_o*theta + C_1;
    surf(x,y,z, 'EdgeColor', 'none')
    hold on

    if (mod(i, Nr) == 0)
        end_theta = start_theta;
        start_theta = start_theta - pi/No;
    end
end

subplot(1,2,2)
ylabel('FE Solution')
start_theta = pi - pi/No;
end_theta = pi;
for i = 1:length(mat(:,1))
    [r, theta] = meshgrid(ri:0.1:ro, linspace(start_theta, end_theta, length(ri:0.1:ro)));
    x = r .* cos(theta);
    y = r .* sin(theta);
    z = 10 * sin(2*theta)/k_th + C_o*theta + C_1;
    surf(x,y,z, 'EdgeColor', 'none')
    hold on

    if (mod(i, Nr) == 0)
        end_theta = start_theta;
        start_theta = start_theta - pi/No;
    end
end


end


