% calls for a new solution to be computed from a sample
% -----
% sol = solution_identify(sample,solver_opt)
% -----
% Input
% -----
% sample = sample object
% solver_opt = options for identifying solution
% ------
% Output
% ------
% sol = solution object
function [sol, solver_opt] = solution_identify(sample,solver_opt)
    sample = folds_identify(sample, solver_opt);
    if ~isfield(solver_opt,'sig')
        solver_opt.sig = 1; % Initialize relative error to 1 for solver, as for using a zero surrogate.
    end
    if  solver_opt.cross_validate % cross-validate
        sol = cross_validate_coefficients(sample,solver_opt);
        solver_opt.sig = min(1,max(sol.sig,sol.err)); % sol.sig is better, but if sol.err larger, indicates higher relative errors. Should always be less than 1
    else % no cross-validation
        sol = validate_coefficients(sample,solver_opt);
    end
end
