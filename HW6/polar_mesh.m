function [coordinates, LM] = polar_mesh(No, Nr, dt, num_nodes, ri, ro, num_elem)

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

LM = zeros(num_elem, 4);
i = 1;
j = i + (Nr + 1);
l = 1; % index for the layer
for o = 1:No
    for r = l:(l + Nr - 1)
        LM(r,:) = [i, j, j + 1, i + 1];
        i = i + 1;
        j = i + (Nr + 1);
    end
    i = i + 1;
    j = i + (Nr + 1);
    l = l + Nr;
end

end