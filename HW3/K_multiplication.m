function [result] = K_multiplication(K, a)
% This function performs multiplication of a matrix K with a vector a

result = zeros(1, length(K(:,1)));
for i = 1:length(K(:,1))
    for j = 1:length(K(1,:))
    result(i) = result(i) + K(i,j) * a(j);
    end
end

result = result';

end

