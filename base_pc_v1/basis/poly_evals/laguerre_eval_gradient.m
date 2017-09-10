% evaluates set of orthonormal Laguerre polynomials
% orthogonal with respect to  gamrnd(a+1,1)
% -----
% f_o = laguerre_eval_gradient(ord,f_i);
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
% fp: rows of function derivative evaluations, column is order
function [f,fp] = laguerre_eval_gradient(ord,f_i,a)
    n = size(f_i,1);
    f = zeros(n,ord);
    fp = zeros(n,ord);
    % These combinations can be used frequently
    a1 = a+1; a2 = a+2; a_1 = a-1;
    switch ord % before 3 term is effective
        case 0 % this case should not be used
            return
        case 1
            f(:,1) = 1+a-f_i;
            fp(:,1) = -1*ones(n,1);
            return
    end
    f(:,1) = 1+a-f_i;
    f(:,2) = 0.5*f_i.^2-a2*f_i+0.5*a2*a1;
    for k = 3:ord % three term recurrence
        f(:,k) = (2*k+a_1-f_i).*f(:,k-1)-(k+a_1)*f(:,k-2);
        f(:,k) = f(:,k)/k;
    end
    fp(:,1) = 1; % forming poly for alpha+1 by sum equality
    for k = 2:ord
        fp(:,k)= fp(:,k-1)+f(:,k-1);
    end
    fp = -fp; % minus sign completes derivative computation
    nc = sqrt(gamma(a1));
    fact = 1;
    for k = 1:ord % normalization
        fact = fact*k;
        c=nc*sqrt(fact/gamma(k+a1));
        f(:,k) = f(:,k)*c;
        fp(:,k) = fp(:,k)*c;
    end
end
