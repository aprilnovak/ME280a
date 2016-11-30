function [var_physical_domain] = PhysicalInterpolation(physical_domain, space_blocks, var_blocks)

var_physical_domain = zeros(1,length(physical_domain));

j = 1;
for i = 1:length(physical_domain)
    
    if (physical_domain(i) <= space_blocks(j))
    else
        j = j + 1;
    end
    
    var_physical_domain(i) = var_blocks(j);
end

end

