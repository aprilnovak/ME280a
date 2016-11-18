function [coordinates] = polar_mesh(No, Nr, dt, num_nodes, ri, ro)

coordinates = zeros(num_nodes, 2);
theta = pi;
c = 1;  % index for coordinate number

for O = 1:(No + 1)

    % mesh in each radial slice
    for r = 1:(Nr + 1)
        coordinates(c,1) = (ri + dt * (r-1)) * cos(theta);
        coordinates(c,2) = -(ri + dt * (r-1)) * sin(theta);
        c = c + 1;
    end

    % increment theta for next slice
    theta = theta + pi / No;

end

end