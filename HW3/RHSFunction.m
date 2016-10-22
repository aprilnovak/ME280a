function [RHS_value] = RHSFunction(gamma, k_freq, space_blocks, E_blocks, node_number)
% function to repeatedly evaluate \circled{A} in the document

loc = space_blocks(node_number);
RHS_value = - loc * (k_freq^3) * cos(gamma * loc) / (gamma ^2) + 2 * k_freq ^ 3 * sin(gamma * loc) / (gamma ^ 3);

end

