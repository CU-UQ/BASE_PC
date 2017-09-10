% adjust basis set by contracting and expanding
% -----
% basis_new = basis_adjust(basis_old, basis_opt, eval_opt)
% -----
% Input
% -----
% basis_old = old basis
% c = coefficients for solution
% tol = tolerance in solution computation
% basis_opt = options associated with basis
% eval_opt = options related to evaluation

% ------
% Output
% ------
% basis_new = adjusted basis object
% basis_opt = new basis_opt values
% sample_rate_adj = flag returns true if sample_rate needs adjustment

function [min_basis, min_basis_opt] = basis_update(basis_old, sample, surrogate, tol, basis_opt, eval_opt, solver_opt, n_iters, max_elems)
    % initialize with initial basis
    min_err = surrogate.err;
    min_basis = basis_old;
    basis = basis_old;
    min_basis_opt = basis_opt;
    old_coeff = basis_opt.expand_coeff; % Lowered each iteration to attempt to identify less expansive expansions.
    old_dim_add = basis_opt.dim_add;
    for k = 1:n_iters        
        % start with contraction amongst minimizing basis
        basis = basis_contract(basis, surrogate.c, tol, basis_opt.remove_perc);
        sample.lhs = basis_eval(basis,sample.rv,eval_opt);
        sample.wlhs = apply_weights(sample.w,sample.lhs);
        surrogate = id_coefficients(sample, solver_opt);        
        % if error reduced over previous basis
        if surrogate.err < min_err
            min_basis = basis;
            min_basis_opt = basis_opt;
            min_err = surrogate.err;
        end
        var_coefs = variance_by_coord(basis,surrogate.c,1);
        var_ord = order_percentage(basis,1);
        ord_indices = max(var_coefs-var_ord,0);
        ord_indices = basis_opt.expand_coeff*ord_indices/max(ord_indices);
        % then expand new basis
        [basis, basis_opt] = basis_coord_expand(basis, basis_opt, eval_opt, ord_indices);
        if basis.n_elems > max_elems % occassionaly the number of elements can be too large.
            % Reduce basis to maximal size
            basis.ord = basis.ord(1:max_elems);
            basis.dim = basis.dim(1:max_elems);
            basis.n_elems = max_elems;
        end
        sample.lhs = basis_eval(basis,sample.rv,eval_opt);
        sample.wlhs = apply_weights(sample.w,sample.lhs);
        surrogate = id_coefficients(sample, solver_opt);
        % check again
        if surrogate.err < min_err
            min_basis = basis;
            min_basis_opt = basis_opt;
            min_err = surrogate.err;
        end
        % reduce tolerance for subsequent iterations
        tol = tol/3;
        % reduce expansion coefficient
        basis_opt.expand_coeff = sqrt(basis_opt.expand_coeff);
        % reduce new dimensions to add per iteration
        basis_opt.dim_add = max(1,ceil(basis_opt.dim_add/1.3));
    end
    min_basis_opt.dim_add = old_dim_add;
    basis_opt.expand_coeff = old_coeff;
end
