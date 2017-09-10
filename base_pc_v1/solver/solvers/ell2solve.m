% Solves least squares, using mldivide
% -----
% x = ell2solve(A,b,inv_lambda)
% -----
% Input
% -----
% A = lhs matrix
% b = rhs vector
% lambda = Tikhonov regularization parameter that may be cross-validated
% ------
% Output
% ------
% x = solution coeffs
function x = ell2solve(A,b,lambda)
if lambda > 0
    nCols = size(A,2);
    A = [A; lambda*eye(nCols)];
    b = [b; zeros(nCols,1)];
end
x = A\b;
end
