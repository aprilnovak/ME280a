function [a_u_condensed, pcg_error] = PCG_better(K, F_u, K_uk, dirichlet_nodes, pcg_error_tol, precondition)

% This function solves the system K * a = R without using the statically-
% condensed matrices.

R = F_u - K_uk * dirichlet_nodes(2,:)';

soln_iter = ones(1, length(K))';

if (precondition == 'precondition')
    % diagonal preconditioner
    T = diag(ones(1,length(K)));
    for i = 1:length(K(1,:))
        T(i,i) = 1 ./ sqrt(K(i,i));
    end

    % transform to the preconditioned regime
    K = transpose(T) * K * T;
    R = transpose(T) * AddDirichlet(R, dirichlet_nodes);
end

%size(R)
R = CutoffDirichlet(R, dirichlet_nodes);
%size(R)
%size(K_multiplication(K, soln_iter, dirichlet_nodes))
r = R - K_multiplication(K, soln_iter, dirichlet_nodes);
z = r;

% compute the initial lambda
lambda = transpose(z) * r / (transpose(z) * K_multiplication(K, AddDirichlet(z, dirichlet_nodes), dirichlet_nodes, non_dn));

% store the previous iteration for convergence estimates
soln_prev = soln_iter;

j = 1;
% perform the first update (only need to add the update term to the
% non-Dirichlet nodes)
for i = 1:length(soln_iter)
    if (find(dirichlet_nodes(1,:) == i))
        soln_iter(i) = soln_prev(i);
    else
        soln_iter(i) = soln_prev(i) + lambda * z(j);
        j = j + 1;
    end
end

j = 1;
% error should only be evaluated with the non-Dirichlet nodes
pcg_error(j) = transpose(CutoffDirichlet(soln_iter - soln_prev, dirichlet_nodes)) * K_multiplication(K, soln_iter - soln_prev, dirichlet_nodes);

num_updates = 1;
while pcg_error(j) > pcg_error_tol
    z_prev = z;
    soln_prev = soln_iter;

    K_soln_prev = K_multiplication(K, soln_prev, dirichlet_nodes);
    K_z_prev = K_multiplication(K, AddDirichlet(z_prev, dirichlet_nodes), dirichlet_nodes);
    
    r = R - K_soln_prev;
    theta = - transpose(r) * K_z_prev / (transpose(z_prev) * K_z_prev);

    z = r + theta * z_prev;
    lambda = transpose(z) * r / (transpose(z) * K_multiplication(K, AddDirichlet(z, dirichlet_nodes), dirichlet_nodes));
    
    j = 1;
    % perform update (only need to add the update term to the non-Dirichlet nodes)
    for i = 1:length(soln_iter)
        if (find(dirichlet_nodes(1,:) == i))
            soln_iter(i) = soln_prev(i);
        else
            soln_iter(i) = soln_prev(i) + lambda * z(j);
            j = j + 1;
        end
    end
    
    K_soln_iter = K_multiplication(K, soln_iter, dirichlet_nodes);
    
    pcg_error(j+1) = transpose(CutoffDirichlet(soln_iter - soln_prev, dirichlet_nodes)) * (K_soln_iter - K_soln_prev);
    num_updates = num_updates + 1;
    j = j + 1;
    
end

if (precondition == 'precondition')
    % transform back from the preconditioned regime
    soln_iter = T * soln_iter;
end

a_u_condensed = CutoffDirichlet(soln_iter, dirichlet_nodes);

sprintf('Number of PCG iterations: %i', num_updates)
end



