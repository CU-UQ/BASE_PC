% validates new basis set from computed candidates
% -----
%[basis_new, basis_opt, sample_rate_adj] = basis_validate(basis_old, c, tol, basis_opt, eval_opt)
% -----
% Input
% -----
% basis_old = old basis
% sample = current sample
% surrogate = surrogate information
% tol = tolerance in solution computation
% basis_opt = options associated with basis
% eval_opt = options related to evaluation

% ------
% Output
% ------
% basis_new = adjusted basis object
% basis_opt = new basis_opt values
% sample_rate_adj = flag returns true if sample_rate needs adjustment

function [basis_new, basis_new_opt] = basis_adjust_with_validation(basis, sample, surrogate, basis_opt, eval_opt)
    if norm(surrogate.c) == 0 % Here, c == 0, so no basis shape information.
        fprintf('Solution is zero vector.\n');
        basis_new = basis;
        basis_new_opt = basis_opt;
        return
    end
    
    % First remove all zero coefficient elements
    if size(surrogate.c,1) == 1 % no contraction from one element
       fprintf('Basis has single function');
       [basis_new, basis_new_opt] = basis_expand(basis,basis_opt,eval_opt); % Expand dand hope fixes problem.
        return
    end
    [c_sorted, indices] = sort(abs(surrogate.c)); % Sort c
    basis.ord = basis.ord(indices); % Reorder these indices
    basis.dim = basis.dim(indices); % Reorder these indices
    ind = find(c_sorted>0,1); % First index
    n_indices = size(indices,1); % total number of indices
    remove_indices = 1:(ind-1); % Remove everything with zero indices
    
    % update basis to remove zero coefficients
    if ~isempty(remove_indices) % Sometimes nothing is removed
        active_indices = ind:basis.n_elems; % We remove first ind-1 elements  
        basis.ord = basis.ord(active_indices);
        basis.dim = basis.dim(active_indices);
        basis.n_elems =  basis.n_elems - ind + 1;
        basis = basis_max_identify(basis);
    end    
    % This is our starting point. From here, we need only track when the
    % basis_max is lowered, inducing a change in basis_new.
    [basis_new, basis_new_opt] = basis_expand(basis, basis_opt, eval_opt); % This is our new expanded basis
    sample.lhs = basis_eval(basis_new,sample.rv,eval_opt); % We don't reset weights as not correction sampling yet
    sample.wlhs = apply_weights(sample.w,sample.lhs); % So we just apply known weights
    basis_exp = basis_new;
    n_dims = size(basis_new.max,2);
    basis_exp_opt = basis_new_opt;
    sol = solution_identify(sample,solver_opt); % Solution in error in this basis
    val_error = sol.err;
    
    % We now remove smallest coefficients until number of strikes is
    % reached or all indices are removed
    strikes = 0;
    cur_max = basis.max;
    cur_max_ord = basis.max_ord;
    while strikes <= basis_opt.strikes && ind < n_indices % 3 strikes typical for baseball analogy Loop over indices
        cur_ord = basis.ord{ind};
        cur_dim = basis.dim{ind};
        basis.ord = basis.ord([1:ind-1 ind+1:basis.n_elems]);
        basis.dim = basis.dim([1:ind-1 ind+1:basis.n_elems]);
        basis.n_elems = basis.n_elems - 1;
        if sum(cur_ord == cur_max(cur_dim)) || sum(cur_ord) == cur_max_ord % Maximums may have changed with this element removed
            basis = basis_max_identify(basis); % Identify new maximum
            if sum(basis.max ~= cur_max) || basis.max_ord ~= cur_max_ord % One will hold if max is changed
                dif_max = basis.max - cur_max; % to adjust expanded basis
                cur_max = basis.max; % Whatever ord has changed
                cur_max_ord = basis.max_ord; % Whatever ord has changed
                dif_max = [dif_max,zeros(1,n_dims- size(dif_max,2))]; %#ok<AGROW> % Extend dif_max to new value
                basis_exp.max = basis_exp.max+dif_max;
                basis_exp.max_ord = max(basis_exp.max);
                n_elems = basis_exp.n_elems;
                removed_indices = [];
                for k = 1: n_elems % Remove all appropriate basis functions from expanded basis
                    cur_ord = basis_exp.ord{k};
                    cur_dim = basis.dim{k};
                    if sum(cur_ord./basis_exp.max(cur_dim)) > 1
                        removed_indices = [removed_indices k]; %#ok<AGROW>
                        basis_exp.ord = basis_exp.ord([1:k-1 k+1:basis_exp.n_elems]);
                        basis_exp.dim = basis_exp.dim([1:k-1 k+1:basis_exp.n_elems]);
                        basis_exp.n_elems = basis_exp.n_elems - 1;
                    end
                end
                if ~isempty(removed_indices) % With removed indices from expanded basis, compute new cv solution
                    indices_to_keep = setdiff(1:n_elems, removed_indices);
                    sample.lhs = sample.lhs(:,indices_to_keep); % As not fixing sampling yet, we just remove columns from lhs and wlhs
                    sample.wlhs = sample.wlhs(:,indices_to_keep);
                    if basis_exp.n_elems <= max_elems % Sometimes we need to contract to get under minimal element size
                        sol = solution_identify(sample,solver_opt); % Solution in error in this basis
                        if (sol.err < val_error) % New minimizing basis
                            basis_new = basis_exp;
                            basis_exp_opt.ord = basis_exp.max;
                            basis_new_opt = basis_exp_opt;
                            val_error = sol_error; % New minimum error to beat
                        else
                            strikes = strikes+1; % Strike because smaller basis didn't reduce error
                        end
                    end
                end
            end
        end
        ind = ind+1; % Prepare to Remove next index
    end
end
