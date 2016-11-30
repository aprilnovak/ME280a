function [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(xe, eta, num_nodes_per_elem, coordinates, LM, elem)

N = zeros(num_nodes_per_elem, 1);
dN = zeros(num_nodes_per_elem, 1);

switch num_nodes_per_elem
    case 4
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
        B = [dN_dxe(1), dN_dxe(2), dN_dxe(3), dN_dxe(4); dN_deta(1), dN_deta(2), dN_deta(3), dN_deta(4)];
    otherwise
        disp('You entered an unsupported number of nodes per element.');
end

% x(xe, eta) and y(xe, eta) transformation to the parametric domain
x_xe_eta = 0.0;
y_xe_eta = 0.0;
dx_dxe = 0.0;
dy_dxe = 0.0;
dx_deta = 0.0;
dy_deta = 0.0;

for i = 1:num_nodes_per_elem
    x_xe_eta = x_xe_eta + coordinates(LM(elem, i), 1) * N(i);
    y_xe_eta = y_xe_eta + coordinates(LM(elem, i), 2) * N(i);
    
    dx_dxe   = dx_dxe   + coordinates(LM(elem, i), 1) * dN_dxe(i);
    dy_dxe   = dy_dxe   + coordinates(LM(elem, i), 2) * dN_dxe(i);
    
    dx_deta = dx_deta   + coordinates(LM(elem, i), 1) * dN_deta(i);
    dy_deta = dy_deta   + coordinates(LM(elem, i), 2) * dN_deta(i);
end