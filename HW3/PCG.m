function [a_u_condensed, a_u_condensed_ge] = PCG(K_uu, F_u, K_uk, dirichlet_nodes, pcg_error_tol)

% This function solves the system K * a = R

K = K_uu;
R = F_u - K_uk * dirichlet_nodes(2,:)';

% pick the initial guess for the solution
soln_iter = ones(1, length(K_uu))';

% compute the initial z from the residual
r = R - K * soln_iter;
z = r;

% compute the initial lambda
lambda = transpose(z) * r / (transpose(z) * K * z);

% store the previous iteration for convergence estimates
soln_prev = soln_iter;

% perform the first update
soln_iter = soln_prev + lambda * z;



pcg_error = 1;
% perform all subsequent updates
while pcg_error > pcg_error_tol
    z_prev = z;
    soln_prev = soln_iter;

    r = R - K * soln_iter;
    theta = - transpose(r) * K * z_prev / (transpose(z_prev) * K * z_prev);

    z = r + theta * z_prev;
    lambda = transpose(z) * r / (transpose(z) * K * z);
    soln_iter = soln_prev + lambda * z;

    pcg_error = (transpose(soln_iter) - transpose(soln_prev)) * K * (soln_iter - soln_prev);
end

% Gaussian elimination method (for comparison)
a_u_condensed_ge = K_uu \ R;

a_u_condensed = soln_iter;

end

