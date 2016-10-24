function [a_u_condensed] = PCG_best(F_u, K_uk, K_cell, LM, num_nodes, dirichlet_nodes, num_elem, num_nodes_per_element, pcg_error_tol, precondition)

% This function solves the system K * a = R

R = F_u - K_uk * dirichlet_nodes(2,:)';

% if (precondition == 'precondition')
%     % diagonal preconditioner
%     T = diag(ones(1, num_nodes));
%     for i = 1:num_nodes
%         T(i,i) = 1 ./ sqrt(K(i,i));
%     end
% 
%     % transform to the preconditioned regime
%     K = transpose(T) * K * T;
%     R = transpose(T) * R;
% end

% pick the initial guess for the solution
soln_iter = ones(1, num_nodes)';

% compute the initial z from the residual
r = R - Mult(K_cell, num_elem, soln_iter, num_nodes, num_nodes_per_element, dirichlet_nodes, LM);
z = r;

% compute the initial lambda
lambda = transpose(z) * r / (transpose(z) * Mult(K_cell, num_elem, AddDirichlet(z, dirichlet_nodes), num_nodes, num_nodes_per_element, dirichlet_nodes, LM));

% store the previous iteration for convergence estimates
soln_prev = soln_iter;

% perform the first update
soln_iter = soln_prev + lambda * AddDirichlet(z, dirichlet_nodes);


j = 1;
pcg_error(j) = transpose(CutoffDirichlet(soln_iter - soln_prev, dirichlet_nodes)) * Mult(K_cell, num_elem, soln_iter - soln_prev, num_nodes, num_nodes_per_element, dirichlet_nodes, LM);
num_updates = 1;

% perform all subsequent updates
while pcg_error(j) > pcg_error_tol
    z_prev = z;
    soln_prev = soln_iter;

    r = R - Mult(K_cell, num_elem, soln_iter, num_nodes, num_nodes_per_element, dirichlet_nodes, LM);
    K_z_prev = Mult(K_cell, num_elem, AddDirichlet(z_prev, dirichlet_nodes), num_nodes, num_nodes_per_element, dirichlet_nodes, LM);
    theta = - transpose(r) * K_z_prev / (transpose(z_prev) * K_z_prev);

    z = r + theta * z_prev;
    lambda = transpose(z) * r / (transpose(z) * Mult(K_cell, num_elem, AddDirichlet(z, dirichlet_nodes), num_nodes, num_nodes_per_element, dirichlet_nodes, LM));
    soln_iter = soln_prev + lambda * AddDirichlet(z, dirichlet_nodes);

    pcg_error(j+1) = transpose(CutoffDirichlet(soln_iter - soln_prev, dirichlet_nodes)) * Mult(K_cell, num_elem, soln_iter - soln_prev, num_nodes, num_nodes_per_element, dirichlet_nodes, LM);
    num_updates = num_updates + 1;
    j = j + 1;
end

if (precondition == 'precondition')
    % transform back from the preconditioned regime
    a_u_condensed = T * a_u_condensed;
end

a_u_condensed = CutoffDirichlet(soln_iter, dirichlet_nodes);

sprintf('Number of PCG iterations: %i', num_updates)



end

