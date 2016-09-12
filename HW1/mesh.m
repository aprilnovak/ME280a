% FE HW 1, April Novak

function [LM, coordinates] = mesh(L, num_elem, shape_order)

% for evenly-spaced nodes, on a 3-D mesh. Each row in coordinates
% corresponds to a node. The node number refers to the row number in the
% coordinates matrix.
coordinates = zeros(num_elem, 3);

% the number of nodes depends on the shape function order
num_nodes = (shape_order - 1) * num_elem + 1;

%in 1-D, the first node starts at (0,0), and the rest are evenly-spaced
for i = 2:num_nodes
   coordinates(i,:) = [coordinates(i - 1, 1) + L/num_elem, 0, 0];
end

% Which nodes correspond to which elements depends on the shape function
% used. Each row in the LM corresponds to one element.
num_nodes_per_element = shape_order;

LM = zeros(num_elem, num_nodes_per_element); 

for i = 1:num_elem
    for j = 1:num_nodes_per_element
        LM(i,j) = num_nodes_per_element * (i - 1) + j;
    end
end

end