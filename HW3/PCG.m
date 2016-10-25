function [a_u_condensed, pcg_error] = PCG(K_uu, F_u, K_uk, dirichlet_nodes, pcg_error_tol, precondition)

K = K_uu;
R = F_u - K_uk * dirichlet_nodes(2,:)';

% pick the initial guess for the solution
soln_iter = ones(1, length(K_uu))';

if (precondition == 'precondition')
    % diagonal preconditioner
    T = diag(ones(1,length(K)));
    for i = 1:length(K(1,:))
        T(i,i) = 1 ./ sqrt(K(i,i));
    end

    % transform to the preconditioned regime
    K = transpose(T) * K * T;
    R = transpose(T) * R;
end

% compute the initial z from the residual
r = R - K * soln_iter;
z = r;

% compute the initial lambda
lambda = transpose(z) * r / (transpose(z) * K * z);

% store the previous iteration for convergence estimates
soln_prev = soln_iter;

% perform the first update
soln_iter = soln_prev + lambda * z;

j = 1;
pcg_error(j) = transpose(soln_iter - soln_prev) * K * (soln_iter - soln_prev);%/(transpose(soln_prev) * K * soln_prev);
num_updates = 1;
% perform all subsequent updates
while pcg_error > pcg_error_tol
    z_prev = z;
    soln_prev = soln_iter;

    r = R - K * soln_iter;
    K_z_prev = K * z_prev;
    theta = - transpose(r) * K_z_prev / (transpose(z_prev) * K_z_prev);

    z = r + theta * z_prev;
    lambda = transpose(z) * r / (transpose(z) * K * z);
    soln_iter = soln_prev + lambda * z;
    
    pcg_error(j+1) = transpose(soln_iter - soln_prev) * K * (soln_iter - soln_prev);% / (transpose(soln_prev) * K * soln_prev);
    num_updates = num_updates + 1;
    j = j + 1;
end

a_u_condensed = soln_iter;

if (precondition == 'precondition')
    % transform back from the preconditioned regime
    a_u_condensed = T * a_u_condensed;
end

sprintf('Number of PCG iterations: %i', num_updates)

end

