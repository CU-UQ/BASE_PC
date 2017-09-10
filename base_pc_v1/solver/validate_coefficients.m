% Computes solution and error estimate without cross-validating tolerance (uses specified tolerance parameter)
% -----
% sol = validate_coefficients(sample,solver_opt)
% -----
% Input
% -----
% sample = sample object
% solver_opt = options for identifying solution
% ------
% Output
% ------
% sol = solution object
function sol = validate_coefficients(sample,solver_opt)
    index = 1:solver_opt.n_workers:sample.n_folds;
    % Main computation
    spmd (solver_opt.n_workers) % Loop over folds to compute error
        f_tol = sol_validate_serial(sample, sample.folds(index+labindex-1,:), solver_opt, solver_opt.sig); % No cross-validation of tolerance
    end
    % Combine the folds
    f_tol = sum(cell2mat(f_tol(:))); % Realizations of ratio of unexplained variance
    sol.c = solver_opt.solver_handle(sample.wlhs, sample.wrhs, solver_opt.sig);
    sol.err = sqrt(f_tol/(solver_opt.n_workers*solver_opt.n_folds_pw));
    sol.sig = solver_opt.sig;
end

function f_tol  = sol_validate_serial(sample, folds, solver_opt, tol)
    n_rows = floor(size(sample.wrhs,1)/sample.n_samps); % For gradient information inclusion
    c_samps = min(floor(solver_opt.comp_samp_perc*sample.n_samps),sample.n_samps-1); % Need at least 1 validation sample.
    v_samps = sample.n_samps - c_samps; % validation samples
    f_tol = 0; % Holding TVA for each tolerance
    for k = 1:solver_opt.n_folds_pw % Loop over number of folds
        this_fold = folds(k,:);
        comp_samps = zeros(n_rows*c_samps,1);% Acquire computation samples
        for j = 1:c_samps % Loop over the appropriate indices
            comp_samps(((j-1)*n_rows+1):(j*n_rows)) = (1+(this_fold(j)-1)*n_rows):(this_fold(j)*n_rows);
        end
        val_samps = zeros(n_rows*v_samps,1); % Acquire validation samples
        for j = 1:v_samps % Loop over the appropriate indices
            val_samps(((j-1)*n_rows+1):(j*n_rows)) = (1+(this_fold(j+c_samps)-1)*n_rows):(this_fold(j+c_samps)*n_rows);
        end
        c_hat = solver_opt.solver_handle(sample.wlhs(comp_samps,:), sample.wrhs(comp_samps), tol); % Compute solution from computation samples
        f_tol = f_tol + (norm(sample.wlhs(val_samps,:)*c_hat-sample.wrhs(val_samps))/norm(sample.wrhs(val_samps)))^2; % This sums estimates of Fraction of Variance Unexplained, computed from validation samples
    end
end
