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

function [basis_new, basis_new_opt, min_solver_opt, min_sample, min_sol, exit_flag] =  basis_validation_anisotropic_total_order(basis, sample, sol, basis_opt, eval_opt, solver_opt)
    original_error = sol.err;
    original_sample = sample;
    original_basis = basis;
    original_basis_opt = basis_opt;
    original_solver_opt = solver_opt;
    original_sol = sol;
    exit_flag = false;
    if ~solver_opt.basis_cross_validate % If false can significantly reduce computational cost
        return_to_cv_value = solver_opt.cross_validate; % Will overwrite this, so need to know to set it back;
        solver_opt.cross_validate = false; % Setting to false to speed computation
    end
    if norm(sol.c) == 0 % Here, c == 0, so no basis shape information. Expand and try.
        fprintf('Solution is zero vector.\n');        
        if basis_opt.order_bounds % Update order bounds if appropriate
            basis_opt = basis_order_bounds(basis_opt, eval_opt);
        end
        [basis, basis_opt] = basis_expand(basis, basis_opt, eval_opt); % This is our new expanded basis
        sample.lhs = basis_eval(basis,sample.rv,eval_opt); % We don't reset weights as not correction sampling yet
        sample.wlhs = apply_weights(sample.w,sample.lhs); % So we just apply known weights
        [sol, solver_opt] = solution_identify(sample,solver_opt); % Solution in error in this basis
        basis_new = basis;
        basis_new_opt = basis_opt;
        min_solver_opt = solver_opt;
        min_sample = sample;
        min_sol = sol;
    end
    
    % First remove all zero coefficient elements
    if size(sol.c,1) == 1 % no contraction from one element
       fprintf('Basis has single function.\n');
       basis_opt.ord = basis.max;
       if basis_opt.order_bounds % Update order bounds if appropriate
           basis_opt = basis_order_bounds(basis_opt, eval_opt);
       end
       [basis, basis_opt] = basis_expand(basis, basis_opt, eval_opt); % This is our new expanded basis
       sample.lhs = basis_eval(basis,sample.rv,eval_opt); % We don't reset weights as not correction sampling yet
       sample.wlhs = apply_weights(sample.w,sample.lhs); % So we just apply known weights
       [sol, solver_opt] = solution_identify(sample,solver_opt); % Solution in error in this basis
       basis_new = basis;
       basis_new_opt = basis_opt;
       min_solver_opt = solver_opt;
       min_sample = sample;
       min_sol = sol;
    end
    [c_sorted, indices] = sort(abs(sol.c)); % Sort c
    basis.ord = basis.ord(indices); % Reorder these indices
    basis.dim = basis.dim(indices); % Reorder these indices
    ind = find(c_sorted>0,1); % First index
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
    basis_opt.ord = basis.max;
    basis_opt.dim  = basis.n_dim;
    if basis_opt.order_bounds % Update order bounds if appropriate
        basis_opt = basis_order_bounds(basis_opt, eval_opt);
    end
    [basis_exp, basis_exp_opt] = basis_expand(basis, basis_opt, eval_opt); % This is our new expanded basis

    sample.lhs = basis_eval(basis_exp,sample.rv,eval_opt); % We don't reset weights as not correction sampling yet
    sample.wlhs = apply_weights(sample.w,sample.lhs); % So we just apply known weights
    val_error = inf;
    if basis_exp.n_elems <= basis_opt.max_elems % If not need to contract until under element limit
        [sol, solver_opt] = solution_identify(sample,solver_opt); % Solution in error in this basis
        val_error = sol.err;
        basis_new = basis_exp;
        basis_new_opt = basis_exp_opt;
        min_solver_opt = solver_opt;
        min_sample = sample;
        min_sol = sol;
    end

    % We now remove smallest coefficients until number of strikes is
    % reached or all indices are removed
    strikes = 0;
    cur_max = basis.max;
    cur_n_dim = basis.n_dim;
    while strikes < basis_opt.validation_strikes && 1 < basis.n_elems % 3 strikes typical for baseball analogy Loop over indices
        cur_ord = basis.ord{1}; % order associated with smallest coefficient
        cur_dim = basis.dim{1}; % dim associated with smallest coefficient
        basis.ord = basis.ord(2:basis.n_elems); % removed from basis
        basis.dim = basis.dim(2:basis.n_elems); % removed from basis
        basis.n_elems = basis.n_elems - 1;
        if sum(cur_ord == cur_max(cur_dim)) % Maximums may have changed with this element removed
            basis = basis_max_identify(basis); % Identify new maximum
            if basis.n_dim < cur_n_dim || sum(basis.max ~= cur_max) % Sometimes dimension is reduced, so evaluate that first or second will be error
                cur_max = basis.max; % Adjust these for future use
                cur_n_dim = basis.n_dim;
                basis_opt.ord = basis.max;
                basis_opt.dim  = basis.n_dim;
                if basis_opt.order_bounds % Update order bounds if appropriate
                    basis_opt = basis_order_bounds(basis_opt, eval_opt);
                end
                basis_exp_opt = basis_expand_opt_identify(basis, basis_opt, eval_opt); % Id new basis opts.
                min_dim = min(basis_exp.n_dim,basis_exp_opt.dim); % Make appropriate changes to the basis, decreasing orders if they are smaller (max elems concerns complicates this slightly)
                basis_exp.max = basis_exp.max(1:min_dim);
                basis_exp.n_dim = min_dim;
                for k = 1: min_dim % Set maximum values
                    basis_exp.max(k) = min(basis_exp_opt.ord(k),basis_exp.max(k)); 
                end
                basis_exp.max_ord = max(basis_exp.max);
                removed_indices = []; 
                for k = 1: basis_exp.n_elems; % Remove all appropriate basis functions from expanded basis
                    cur_ord = basis_exp.ord{k};
                    cur_dim = basis_exp.dim{k};
                    if ~isempty(cur_dim) % We end up testing the constant function which messes up the next test
                        if (max(cur_dim) > basis_exp.n_dim) || (sum(cur_ord./basis_exp.max(cur_dim)) > 1) % Again we test the condition first because 2nd will be error
                          removed_indices = [removed_indices k]; %#ok<AGROW> Will remove that basis function
                        end
                    end
                end
                if ~isempty(removed_indices) % With removed indices from expanded basis, compute new cv solution. If max elems was hit, it's possible removed indices is blank.
                    kept_indices = setdiff(1:basis_exp.n_elems, removed_indices);
                    basis_exp.ord = basis_exp.ord(kept_indices);
                    basis_exp.dim = basis_exp.dim(kept_indices);
                    basis_exp.n_elems = size(basis_exp.ord,1);
                    sample.lhs = sample.lhs(:,kept_indices); % As not fixing sampling yet, we just remove columns from lhs and wlhs
                    sample.wlhs = sample.wlhs(:,kept_indices);
                    if basis_exp.n_elems <= basis_opt.max_elems % Sometimes we need to contract to get under minimal element size
                        [sol, solver_opt] = solution_identify(sample,solver_opt); % Solution in error in this basis
                        if (sol.err < val_error) % New minimizing basis and options. Only test if valid n_elems size
                            basis_new = basis_exp;
                            basis_new_opt = basis_exp_opt;
                            val_error = sol.err; % New minimum error to beat
                            min_solver_opt = solver_opt;
                            min_sample = sample;
                            min_sol = sol;
                            strikes = 0; % Reset strikes with identification of new minimum.
                        else
                            strikes = strikes+1; % Strike count if smaller basis didn't reduce error, break if too many strikes
                        end
                    end
                end
            end
        end
    end
    if ~solver_opt.basis_cross_validate % If false can significantly reduce computational cost
        min_solver_opt.cross_validate = return_to_cv_value;
    end
    if val_error >= original_error % Possible that all of these bases were bad
        min_sol = original_sol;
        min_sample = original_sample;
        basis_new = original_basis;
        basis_new_opt = original_basis_opt;
        min_solver_opt = original_solver_opt;
        exit_flag = true;
    end
end
