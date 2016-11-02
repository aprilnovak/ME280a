% Mesh generator, ME 280a HW 5

Nt = 2;         % number of layers
No = 1;         % number of elements in the theta direction
Nc = 8;         % number of elements in circumferential direction
num_nodes_per_elem = 4;     % linear elements

R = 1;          % radius of each arch
r = 0.2;        % radius of inner hole
t = 0.1;        % thickness of the tube wall

layer_thickness = t / (Nt);           % thickness of each ring
angle = (2*pi) / Nc;                % angle in horizontal plane


% in a particular horizontal (relative to tube axis) plane:

% global node numbering begins on the inner surface, and moves clockwise
% until reaching the outer surface

% each row represents one coordinate of a global node
coordinates = zeros(Nt * Nc, 3);


theta = 0;
dt = 0;          

k = 1; % index for coordinate row
for j = 1:(Nt + 1) % create all layers of rings
    
    for i = 1:Nc % create a single ring

        % x-coordinate
        coordinates(k, 1) = (r + dt) * cos(theta);

        % y-coordinate
        coordinates(k, 2) = (r + dt) * sin(theta);

        % z-coordinate
        coordinates(k,3) = 0;
        
        k = k + 1;

        theta = theta + angle;
    end
    
    dt = dt + layer_thickness; % reset radius
    theta = 0; % reset angle
end

x = coordinates(:,1);
y = coordinates(:,2);
z = coordinates(:,3);

scatter3(x,y,z)
