function [a_u_condensed] = PCG_best(F_u, K_uk, K_cell, LM, num_nodes, dirichlet_nodes, num_elem, num_nodes_per_element, pcg_error_tol, precondition)

% figure out how to eliminate the CutoffDirichlet and AddDirichlet
% functions

% This function solves the system K * a = R

R = F_u - K_uk * dirichlet_nodes(2,:)';
R = AddDirichlet(R, dirichlet_nodes);

% if (precondition == 'precondition')
%     % diagonal preconditioner
%     T = diag(ones(1, num_nodes));
%     for i = 1:length(K(1,:))
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
%lambda = transpose(CutoffDirichlet(z, dirichlet_nodes)) * CutoffDirichlet(r, dirichlet_nodes) / (transpose(CutoffDirichlet(z, dirichlet_nodes)) * CutoffDirichlet(Mult(K_cell, num_elem, z, num_nodes, num_nodes_per_element, dirichlet_nodes, LM), dirichlet_nodes));
lambda = transpose(z(2:(end-1))) * r(2:(end-1)) / (transpose(z(2:(end-1))) * CutoffDirichlet(Mult(K_cell, num_elem, z, num_nodes, num_nodes_per_element, dirichlet_nodes, LM), dirichlet_nodes)); 

% store the previous iteration for convergence estimates
soln_prev = soln_iter;

% perform the first update
soln_iter = soln_prev + lambda * z;


j = 1;
% only calculate error using the non-Dirichlet nodes
%disp('New method:')
%transpose(CutoffDirichlet(soln_iter - soln_prev, dirichlet_nodes))
%CutoffDirichlet(Mult(K_cell, num_elem, soln_iter - soln_prev, num_nodes, num_nodes_per_element, dirichlet_nodes, LM), dirichlet_nodes)
K_soln_iter = Mult(K_cell, num_elem, soln_iter - soln_prev, num_nodes, num_nodes_per_element, dirichlet_nodes, LM);
pcg_error(j) = transpose(soln_iter(2:(end-1)) - soln_prev(2:(end-1))) * K_soln_iter(2:(end-1));
num_updates = 1;

% perform all subsequent updates
while pcg_error(j) > pcg_error_tol
    z_prev = z;
    soln_prev = soln_iter;

    r = R - Mult(K_cell, num_elem, soln_iter, num_nodes, num_nodes_per_element, dirichlet_nodes, LM);
    K_z_prev = Mult(K_cell, num_elem, z_prev, num_nodes, num_nodes_per_element, dirichlet_nodes, LM);
    theta = - transpose(r(2:(end-1))) * K_z_prev(2:(end-1)) / (transpose(z_prev(2:(end-1))) * K_z_prev(2:(end-1)));

    z = r + theta * z_prev;
    K_z = Mult(K_cell, num_elem, z, num_nodes, num_nodes_per_element, dirichlet_nodes, LM);
    lambda = transpose(z(2:(end-1))) * r(2:(end-1)) / (transpose(z(2:(end-1))) * K_z(2:(end-1)));
    soln_iter = soln_prev + lambda * z;

    K_soln_iter_soln_prev = Mult(K_cell, num_elem, soln_iter - soln_prev, num_nodes, num_nodes_per_element, dirichlet_nodes, LM);
    pcg_error(j+1) = transpose(soln_iter(2:(end-1)) - soln_prev(2:(end-1))) * K_soln_iter_soln_prev(2:(end-1));
    num_updates = num_updates + 1;
    j = j + 1;
end

% if (precondition == 'precondition')
%     % transform back from the preconditioned regime
%     a_u_condensed = T * a_u_condensed;
% end

a_u_condensed = CutoffDirichlet(soln_iter, dirichlet_nodes);

sprintf('Number of PCG iterations: %i', num_updates)

end

