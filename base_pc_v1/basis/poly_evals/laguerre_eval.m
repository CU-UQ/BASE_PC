% evaluates set of orthonormal Laguerre polynomials
% orthogonal with respect to  gamrnd(a+1,1)
% -----
% f_o = laguerre_eval(ord,f_i);
% -----
% Input
% -----
% ord: maximum order of polynomial
% f_i: column vector of points being evaluated
% a: Laguerre alpha
% ------
% Output
% ------
% f: rows of function evaluations, column corresponds to order
function f = laguerre_eval(ord,f_i,a)
    % These combinations can be used frequently
    a1 = a+1; a2 = a+2; a_1 = a-1;
    f = zeros(size(f_i,1),ord);
    switch ord % before 3 term is effective
        case 0 % this case should not be used
            return
        case 1
            f(:,1) = 1+a-f_i;
            return
    end
    f(:,1) = 1+a-f_i;
    f(:,2) = 0.5*f_i.^2-a2*f_i+0.5*a2*a1;
    for k = 3:ord % three term recurrence
        f(:,k) = (2*k+a_1-f_i).*f(:,k-1)-(k+a_1)*f(:,k-2);
        f(:,k) = f(:,k)/k;
    end
    nc = sqrt(gamma(a1));
    fact = 1;
    for k = 1:ord % normalization
        fact = fact*k;
        f(:,k) = f(:,k)*nc*sqrt(fact/gamma(k+a1));
    end
end
