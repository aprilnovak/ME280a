function [result] = NodalInterpolation(vector, coordinates, physical_domain)
% This function interpolates a vector in the physical domain to a vector
% defined only at nodal values.

j = 1;
k = 1;
for i = 1:length(physical_domain)
    if (abs(physical_domain(i) - coordinates(j)) < 1e-10)
        result(k) = vector(i);
        k = k + 1;
        j = j + 1;
    else
    end
end

end

