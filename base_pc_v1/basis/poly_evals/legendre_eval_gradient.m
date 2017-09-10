% evaluates set of orthonormal Legendre polynomials
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
% fp: rows of function derivative evaluations, column is order
function [f,fp] = legendre_eval_gradient(ord,f_i)
    n = size(f_i,1);
    f = zeros(n,ord);
    fp = zeros(n,ord);
    switch ord
        case 0
            return
        case 1               
            fp(:,1) = ones(n,1)*sqrt(3);
            f(:,1) = f_i*sqrt(3);
            return
        otherwise
            fp(:,1) = ones(n,1);
            fp(:,2) = 3*f_i;
            f(:,1) = f_i;
            f(:,2) = (1.5)*f_i.*f(:,1) - .5;
            for k = 3:ord % recurrence
                f(:,k) = (2-1/k)*f_i.*f(:,k-1) - (1-1/k)* f(:,k-2);
            end
    end
    for k = 3:ord
        fp(:,k) = k*f(:,k-1)+ f_i.*fp(:,k-1);
    end
    for k = 1:ord % normalization
        f(:,k) = f(:,k)*sqrt(2*k+1);
        fp(:,k) = fp(:,k)*sqrt(2*k+1);
    end
end
