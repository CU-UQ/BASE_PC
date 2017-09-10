% Evaluates set of orthonormal Legendre polynomials
% orthogonal with respect to 2*rand-1
% -----
% f_o = legendre_eval(ord,f_i);
% -----
% Input
% -----
% ord: maximum order of polynomial
% f_i: column vector of points being evaluated
% ------
% Output
% ------
% f: rows of function evaluations, column corresponds to order
function f = legendre_eval(ord,f_i)
    f = zeros(size(f_i,1),ord);
    switch ord % before 3 term is effective
        case 0 % this case should not be used
            return
        case 1
            f(:,1) = f_i*sqrt(3);
            return
    end
    f(:,1) = f_i;
    f(:,2) = (3*f_i.^2-1)/2;
    for i = 3:ord % three term recurrence
        t = 1/i;
        f(:,i) = (2-t)*f_i.*f(:,i-1) - (1-t)* f(:,i-2);
    end
    
    for i=1:ord % normalization
        f(:,i) = f(:,i)*sqrt(2*i+1);
    end
end
