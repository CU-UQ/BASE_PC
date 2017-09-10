% evaluates set of orthonormal Hermite polynomials (probabilist)
% orthogonal with respect to randn
% -----
% f_o = hermite_eval(ord,f_i);
% -----
% Input
% -----
% ord: maximum order of polynomial
% f_i: column vector of points being evaluated
% ------
% Output
% ------
% f: rows of function evaluations, column corresponds to order
function f = hermite_eval(ord,f_i)
    f = zeros(size(f_i,1),ord);
    switch ord % before 3 term is effective
        case 0 % this case should not be used
            return
        case 1
            f(:,1) = f_i;
            return
    end
    f(:,1) = f_i;
    f(:,2) = f_i.*f(:,1)-1;
    for k = 3:ord % three term recurrence
        f(:,k) = f_i.*f(:,k-1)-(k-1)*f(:,k-2);
    end
    fact = 1;
    for k = 2:ord % normalization
        fact = fact*k;
        f(:,k) = f(:,k)/sqrt(fact);
    end
end
