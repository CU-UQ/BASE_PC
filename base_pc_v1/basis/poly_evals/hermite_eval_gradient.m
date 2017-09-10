% evaluates set of orthonormal hermite polynomials (probabilist) with derivatives
% orthogonal with respect to randn
% -----
% [f, fp] = hermite_eval_gradient(ord,f_i)
% ----- 
% Input
% -----
% ord: maximum order of polynomial
% f_i: column vector of points being evaluated

% Output
% ------
% f : rows of function evaluations, column corresponds to order
% fp: rows of function derivative evaluations, column is order
function [f, fp] = hermite_eval_gradient(ord,f_i)
    n = size(f_i,1);
    f = zeros(n,ord);
    fp = zeros(n,ord);
    switch ord
        case 0
            return
        case 1
            f(:,1) = f_i;
        otherwise
            f(:,1) = f_i;
            f(:,2) = f_i.*f(:,1)-1;
            for k = 2:ord-1
                f(:,k+1) = f_i.*f(:,k)-k*f(:,k-1);
            end
    end
    %
    fac = 1;
    fp(:,1) = ones(n,1);
    for k = 2:ord
        fac = fac*k;
        f(:,k) = f(:,k)/sqrt(fac);        
        fp(:,k) = sqrt(k)*f(:,k-1);
    end
end
