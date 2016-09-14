% generate x(xe), the shape functions (in the master domain), and the Jacobian

% N - shape functions in the master domain
% dN - derivative of the shape functions with respect to xe
function [N, dN, x_xe, dx_dxe] = shapefunctions(evaluation_point, shape_order, coordinates, LM, elem)

% shape functions and their derivatives
N = zeros(shape_order, 1);
dN = zeros(shape_order, 1);

switch shape_order
    case 2
        N(1) = (1 - evaluation_point) / 2;
        N(2) = (1 + evaluation_point) / 2;
        dN(1) = - 1/2;
        dN(2) = 1/2;
    otherwise
        disp('You entered an unsupported shape function order.');
end

% x(xe) transformation to the parametric domain
x_xe = 0.0;
dx_dxe = 0.0;
for i = 1:shape_order
    x_xe = x_xe + coordinates(LM(elem, i)) * N(i);
    dx_dxe = dx_dxe + coordinates(LM(elem,i)) * dN(i);
end