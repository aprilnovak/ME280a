function [var_elem, right_endpoint_index, right_endpoint_coordinate] = ElementInterpolation(coordinates, num_elem, num_nodes_per_element, space_blocks, var_blocks)

% This function interpolates a block-structued variable into the physical
% domain. By default, values at the ends of boundaries are put into the
% previous domain.

    s = 1;
    var_elem = zeros(1, length(num_elem));
    right_endpoint_index = zeros(1, length(num_elem));
    right_endpoint_coordinate = zeros(1, length(num_elem));
    
    for elem = 1:num_elem
        right_endpoint_index(elem) = elem * (num_nodes_per_element - 1);
        right_endpoint_coordinate(elem) = coordinates(right_endpoint_index(elem) + 1, 1);

        if (abs(right_endpoint_coordinate(elem) - space_blocks(s)) < 1e-10) || (right_endpoint_coordinate(elem) <= space_blocks(s))
            
        else
            s = s + 1;
        end
        var_elem(elem) = var_blocks(s);
    end


end

