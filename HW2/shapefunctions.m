% N         : shape functions in the master domain
% dN        : derivative of the shape functions with respect to xe
% x_xe      : x as a function of xe
% dx_dxe    : derivative of x with respect to xe
function [N, dN, x_xe, dx_dxe] = shapefunctions(xe, shape_order, coordinates, LM, elem)

% shape functions and their derivatives
N = zeros(shape_order, 1);
dN = zeros(shape_order, 1);

switch shape_order
    case 2
        N(1) = (1 - xe) ./ 2;
        N(2) = (1 + xe) ./ 2;
        dN(1) = - 1/2;
        dN(2) = 1/2;
    case 3
        N(1) = xe .* (xe - 1) ./ 2;
        N(2) = - (xe - 1) .* (1 + xe);
        N(3) = xe .* (1 + xe) ./ 2;
        dN(1) = xe - 1/2;
        dN(2) = -2 .* xe;
        dN(3) = 1/2 + xe;
    case 4
        N(1) = (9/16) * (1 - xe) * (1/3 - xe) * (-1/3 - xe);
        N(2) = (-27/16) * (1 - xe) * (1/3 - xe) * (-1 - xe);
        N(3) = (27/16) * (1 - xe) * (-1/3 - xe) * (-1 - xe);
        N(4) = (-9/16) * (1/3 - xe) * (-1/3 - xe) * (-1 - xe);
        dN(1) = (9/16) * (-3 * xe.^ 2 + 2 * xe + 1/9);
        dN(2) = (-27/16) * (-3 * xe .^ 2 + 2 .* xe ./ 3 + 1);
        dN(3) = (27/16) * (-3 * xe .^ 2 - 2.* xe ./ 3 + 1);
        dN(4) = (-9/16) * (-3 * xe ^ 2 - 2 * xe + 1/9);
    otherwise
        disp('You entered an unsupported shape function order.');
end

% check that the sum of the shape functions adds up to 1
sum = 0;
for j = 1:shape_order
    sum = sum + N(j);
end

if (abs(sum - 1.0) > 1e-10)
    disp('Sum of the shape functions does not add up to 1.');
end

% x(xe) transformation to the parametric domain
x_xe = 0.0;
dx_dxe = 0.0;
for i = 1:shape_order
    x_xe = x_xe + coordinates(LM(elem, i)) * N(i);
    dx_dxe = dx_dxe + coordinates(LM(elem,i)) * dN(i);
end