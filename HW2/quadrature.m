function [wt, qp] = quadrature(shape_order)

switch shape_order
    case 2
        wt = [1.0, 1.0];
        qp = [-sqrt(1/3), sqrt(1/3)];
    case 3
        wt = [5/9, 8/9, 5/9];
        qp = [-sqrt(3/5), 0, sqrt(3/5)];
    otherwise
        disp('You entered an unsupported shape function order for the quadrature rule.');
end