% Mesh generator, ME 280a HW 5

Nt = 1;         % number of layers
No = 1;         % number of elements in the theta direction
Nc = 8;         % number of elements in circumferential direction
num_nodes_per_elem = 4;     % linear elements

R = 1;          % radius of each arch
r = 0.2;        % radius of inner hole

% in a particular horizontal (relative to tube axis) plane:

% global node numbering begins on the inner surface, and moves clockwise
% until reaching the outer surface

% each row represents one coordinate of a global node
coordinates = zeros(Nt * Nc * No * num_nodes_per_elem, 3);

angle = (2*pi) / Nc;
theta = 0;
for i = 1:Nc % create a layer
    
    % x-coordinate
    coordinates(i, 1) = r * cos(theta);
    theta = theta + angle;
end




