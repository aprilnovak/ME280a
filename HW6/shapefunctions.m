% N         : shape functions in the master domain
% dN        : derivative of the shape functions with respect to xe
% x_xe      : x as a function of xe
% dx_dxe    : derivative of x with respect to xe
function [N, dN_dxe, dN_deta, x_xe, dx_dxe] = shapefunctions(xe, eta, shape_order, coordinates, LM, elem)

% shape functions and their derivatives (for linear elements)
N = zeros(shape_order * 2, 1);
dN = zeros(shape_order * 2, 1);

switch shape_order
    case 2
        N(1) = 0.25 * (1 - xe) * (1 - eta);
        N(2) = 0.25 * (1 + xe) * (1 - eta);
        N(3) = 0.25 * (1 + xe) * (1 + eta);
        N(4) = 0.25 * (1 - xe) * (1 + eta);
        dN_dxe(1) = 0.25 * -(1 - eta);
        dN_dxe(2) = 0.25 * (1 - eta);
        dN_dxe(3) = 0.25 * (1 + eta);
        dN_dxe(4) = 0.25 * -(1 + eta);
        dN_deta(1) = 0.25 * (1 - xe) * -1;
        dN_deta(2) = 0.25 * (1 + xe) * -1;
        dN_deta(3) = 0.25 * (1 + xe) * 1;
        dN_deta(4) = 0.25 * (1 - xe) * 1;
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