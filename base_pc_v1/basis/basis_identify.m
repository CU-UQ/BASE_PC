% updates basis object
% -----
% basis = basis_identify(basis_opt, eval_opt)
% -----
% Input
% -----
% basis = basis object
% sample = sample object
% sol = solution object
% basis_opt = options to identify basis
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% min_basis = basis object
% min_basis_opt = basis options
% min_solver_opt = solver options
% min_sample = sample object adjusted for new basis
% min_sol = solution object adjusted for new basis
function [min_basis, min_basis_opt, min_solver_opt, min_sample, min_sol] = basis_identify(basis, sample, sol, basis_opt, eval_opt, solver_opt) % can be modified for different bases
    val_err = inf;
    min_basis = basis;
    min_basis_opt = basis_opt;
    min_solver_opt = solver_opt;
    min_sample = sample;
    min_sol = sol;
    for k = 1:basis_opt.validation_iters
        [basis1, basis_opt, solver_opt, sample, sol, exit_flag] = basis_opt.basis_identify_handle(basis, sample, sol, basis_opt, eval_opt, solver_opt);
        if  sol.err < val_err
            val_err = sol.err;
            min_basis = basis1;
            min_basis_opt = basis_opt;
            min_solver_opt = solver_opt;
            min_sample = sample;
            min_sol = sol;
        end
        fprintf('Basis Adaptation Iteration: %i Complete. Validation Error: %f Number of Basis Elements: %i \n', k, val_err, min_basis.n_elems)
        if exit_flag % original basis was selected, future iterations unneeded
            return
        else
            basis = basis1;
        end
    end
    fprintf('Basis Adaptation Complete.\n')
end
