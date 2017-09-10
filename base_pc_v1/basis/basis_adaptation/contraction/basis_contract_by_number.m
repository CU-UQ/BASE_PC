% contracts basis set by specific number
% -----
% basis = basis_contract_by_number(basis, c, number)
% -----
% Input
% -----
% basis = basis object
% c = solution coefficients for corresponding basis
% n_keep = approximate number of basis functions to keep
% ------
% Output
% ------
% basis = contracted basis object

function basis = basis_contract_by_number(basis, c, n_keep)
    n_rem = basis.n_elems - n_keep; % Number of elements to remove
    % First remove all zero coefficient elements
    x = find(c == 0); % set of indices to be removed
    c = abs(c);
    % Second remove smallest coefficients up to tol/3 threshold.
    t_rem = size(x,1);
    while t_rem < n_rem
        i = min(c(c>0));
        if isempty(i)
            break
        end
        j = find(c==i,1);
        t_rem = t_rem + 1;
        if t_rem < n_rem % Err on the side of an extra basis function
            x = [x; j]; %#ok<AGROW>
            c(j) = 0;
        end
    end

    % update basis appropriately
    if ~isempty(x) % Sometimes nothing is removed
        x = setdiff(1:basis.n_elems,x); % We tracked what we want to remove.    
        basis.ord = basis.ord(x);
        basis.dim = basis.dim(x);
        basis.n_elems = size(basis.ord,1);
        basis = basis_max_identify(basis);
    end
end
