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

function [basis, basis_opt, solver_opt, sol] = basis_adjust(basis, sample, sol, basis_opt, eval_opt, solver_opt, max_elems);
    if norm(sol.c) == 0 % Here, c == 0, so no basis shape information.
        fprintf('Solution is zero vector.\n');
        basis_new = basis_old;
        return
    end    
    basis_new = basis_contract(basis_old, sol, sol.err*basis_opt.remove_perc); % Contract basis based on validated error
    if basis_opt.order_bounds
        basis_opt = basis_order_bounds(basis_opt, eval_opt);
    end
    [basis_new, basis_opt] = basis_expand(basis_new, basis_opt, eval_opt);
    if basis.n_elems > max_elems
        basis.ord = basis.ord(1:max_elems);
        basis.dim = basis.dim(1:max_elems);
        basis = basis_max_identify(basis); % Identify new maximum
    end
end
