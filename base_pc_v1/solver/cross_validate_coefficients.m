% Computes solution with cross-validated error
% -----
% sol = cross_validate_coefficients(sample,solver_opt)
% -----
% Input
% -----
% sample = sample object
% solver_opt = options for identifying solution
% ------
% Output
% ------
% sol = solution object
function sol = cross_validate_coefficients(sample,solver_opt)
    index = 1:solver_opt.n_workers:sample.n_folds; % Splits folds appropriately
    % Main computation
    min_log_tol = min(-2,solver_opt.log_tol_exp*log(solver_opt.sig+10^(-10))); % minimum log of tolerance considered here is always a bit lower than desired tolerance
    max_log_tol = min(0, log((solver_opt.log_tol_exp+1)*solver_opt.sig+10^(-10))); % maximum log of tolerance considered. solver_opt.sig is approximately of order of tolerance in ideal world. This gives leeway
    inc = (max_log_tol-min_log_tol)/(solver_opt.cv_mesh_size-2);
    t_tols = [0 exp(min_log_tol:inc:max_log_tol)]'; % Test Tolerances. 
    spmd (solver_opt.n_workers)
        f_tols = sol_cv_serial(sample, sample.folds(index+labindex-1,:), solver_opt, t_tols);
    end
    % Combine the folds
    f_tols = cell2mat(f_tols(:)); % Realizations of ratio of unexplained variance
    n_tols = size(t_tols,1);
    f = zeros(n_tols,1);
    index = 1:n_tols;
    for k = 1:solver_opt.n_workers
        f = f+f_tols(index); % Sums all relevant realized errors over sample.n_folds (Fraction of variance unexplained)
        index = index+n_tols;
    end

    % Best_tol minimizes f over the folds
    [best_f, best_i] = min(f);    
    best_sig = t_tols(best_i);
    
    sol.err = sqrt(best_f/(solver_opt.n_folds_pw*solver_opt.n_workers)); % This is the estimate of RRMSE (Relative Root Mean Square Error)
    sol.sig = best_sig;
    sol.c = solver_opt.solver_handle(sample.wlhs, sample.wrhs, best_sig);

end

function f_tols  = sol_cv_serial(sample, folds, solver_opt, t_tols)
    n_rows = floor(size(sample.wrhs,1)/sample.n_samps); % For gradient information inclusion
    c_samps = min(floor(solver_opt.comp_samp_perc*sample.n_samps),sample.n_samps-1); % Need at least 1 validation sample.
    v_samps = sample.n_samps - c_samps; % validation samples
    n_tols = size(t_tols,1); % Number of tolerances considered
    f_tols = zeros(n_tols,1); % Holding TVA for each tolerance
    for kk = 1:n_tols % Loop over tolerances
        tol = t_tols(kk);
        for k = 1:solver_opt.n_folds_pw % Loop over number of folds
            this_fold = folds(k,:);
            if n_rows > 1
                comp_samps = zeros(n_rows*c_samps,1);% Acquire computation samples
                for j = 1:c_samps % Loop over the appropriate indices
                    comp_samps(((j-1)*n_rows+1):(j*n_rows)) = (1+(this_fold(j)-1)*n_rows):(this_fold(j)*n_rows);
                end
                val_samps = zeros(n_rows*v_samps,1); % Acquire validation samples
                for j = 1:v_samps % Loop over the appropriate indices
                    val_samps(((j-1)*n_rows+1):(j*n_rows)) = (1+(this_fold(j+c_samps)-1)*n_rows):(this_fold(j+c_samps)*n_rows);
                end
            else
                comp_samps = this_fold(1:c_samps); % Conceptually much quicker when n_rows = 1
                val_samps = this_fold(c_samps+1:sample.n_samps);  % Conceptually much quicker when n_rows = 1
            end
            c_hat = solver_opt.solver_handle(sample.wlhs(comp_samps,:), sample.wrhs(comp_samps), tol); % Compute solution vfom computation samples
            f_tols(kk) = f_tols(kk) + (norm(sample.wlhs(val_samps,:)*c_hat-sample.wrhs(val_samps))/norm(sample.wrhs(val_samps)))^2; % This sums estimates of Fraction of Variance Unexplained, computed from validation samples
        end
    end
end
