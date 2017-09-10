% orthogonal matching pursuit
% -----
% x = omp(A,b,tol)
% -----
% Input
% -----
% A = lhs matrix
% b = rhs vector
% tol = regularization parameter that may be cross-validated
% ------
% Output
% ------
% x = coefficient vector
% x = omp(A, b, tol)
% A,x, b = from Ax=b
% tol = For stopping OMP and matrix operations.

% This code modifies SolveOMP.m as part of SparseLAB https://sparselab.stanford.edu/ having the following copyrights
% -----------------------------------------------------
% Copyright (c) 2006. Victoria Stodden and David Donoho
% Copyright (c) 2006. Yaakov Tsaig
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
function x = omp(A,b,tol)
    % Initialize
    eps = 1e-12; % Effective zero is needed due to some numerical instabilities
    tol = max(tol,eps); % if tol is too low, increase
    numRows = size(A,1);
    numCols = size(A,2);
    w = zeros(numCols,1); % weights make correlations most comparable
    for k = 1:numCols
        w(k) = 1/norm(A(:,k));
        A(:,k) = A(:,k)*w(k);
    end  
    x = zeros(numCols,1);
    actSet = []; % active set
    L = []; % For Cholesky decomposition
    nb = norm(b);
    if(nb <= eps) % RHS is numerically zero
        return
    end
    b = b/nb; % Normalize rhs (will renormalize solution later)
    r = b; % Residual to rhs.    
    corrA = A'*b; % Correlation with RHS
    converged = false;
    
    % Main Loop
    iCols=0;
    while(~converged)
        iCols=iCols+1; % Size of ActiveSet        
        corrSense = A'*r;
        corrSense(actSet) = zeros(size(actSet,2),1);
        [C, i] = max(abs(corrSense));
        if C < eps
            x = w.*x*nb;
            return
        end
        [L, actSet, flag] = updateChol(L, A, actSet, i);
        if flag % Gramian has numerical instability
            x = zeros(numCols,1);
            x(actSet) = A(:,actSet)\b; % May not be ok?
	    x = w.*x*nb;
            return
        end        
        % Calculate New Residual
        x = zeros(numCols,1);
        x(actSet) = qcs(L,corrA(actSet));
        % Compute new residual
        r = b - A(:,actSet) * x(actSet);
        % Check Convergence
        if(iCols==numRows-1||norm(r)<tol)
            converged=true;
        end
    end
    x = w.*x*nb;
end
% L = Upper Triangular Cholesky matrix
% LL^T = A^TA. Restricted to columns in active Set
% LL^Tx = A^TAx = z.
function z = qcs(L, z)
    o.LT = true; 
    ot.LT = true; ot.TRANSA = true;
    [w,~] = linsolve(L,z,o);
    [z,~] = linsolve(L,w,ot);
end
% L = existing Lower Triangular Cholesky matrix
% index = index to be added
function [L, activeSet, flag] = updateChol(L, A, activeSet, index) % Updates Cholesky factorization
    flag = 0;
    o.LT = true;
    nCols = size(activeSet,2);
    if(nCols==0)
        L = norm(A(:,index));
    else
        [p,~] = linsolve(L,A(:,activeSet)'*A(:,index),o);
        q = 1-sum(p.^2);
        if q <= 0 % Can occur due to numerical issues
            flag = 1;
            return
        end
        L = [L zeros(nCols,1); p' sqrt(q)];
    end
    activeSet = [activeSet index];
end
