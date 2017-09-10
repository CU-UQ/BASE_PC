% contracts basis set based on solution coefficients
% -----
% basis = basis_contract(basis, c, tol, remove_perc)
% -----
% Input
% -----
% basis = basis object
% c = solution coefficients for corresponding basis
% tol = solution tolerance
% remove_perc = percentage of tolerance for discarding basis functions
% ------
% Output
% ------
% basis = contracted basis object

function basis = basis_contract(basis, c, remove_perc)
    % First remove all zero coefficient elements
    if size(c,1) == 1 % no contraction from one element
        return
    end
    x = find(c == 0);
    c = abs(c);
    % Second remove smallest coefficients up to tol/3 threshold.
    i_sum = 0;
    c = c/norm(c,2); % Normalize so that entries correspond to fraction of variance explained roughly
    t_var = remove_perc; %Err on the side of taking extra basis functions remove coefficients up to removal_perc of tolerances
    while i_sum < t_var
        i = min(c(c>0));
        if isempty(i)
            break
        end
        j = find(c==i,1);
        i_sum = i_sum+i^2; % Sum of variance kept to throw away
        if i_sum < t_var % Err on the side of taking extra basis functions
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
        basis_max_identify(basis);
    end
end
