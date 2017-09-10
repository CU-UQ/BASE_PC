% generates default options from few problem/computer dependent options
% -----
% [basis_opt, eval_opt, sample_opt, solver_opt] = default_parameters(pool_data,dim,eval_type,qoi_handle)
% -----
% Input
% -----
% pool_data = grabs pool_data.NumWorkers for use
% dim = dimension of problem (scalar)
% eval_type = type of random variable (character)
% qoi_handle = handle object for evaluating qoi (@zero_eval is a good choice if QoI not of interest or not called from MATLAB
% ------
% Output
% ------
% basis_opt = options associated with basis
% eval_opt = options associated with evaluation
% sample_opt = options associated with sample
% solver_opt = options associated with solver
function [basis_opt, eval_opt, sample_opt, solver_opt] = default_parameters(pool_data,dim,eval_type,qoi_handle)
    % Slated for removal
    eval_opt.parallel_qoi_eval = false; % Rarely useful for these implementations
    solver_opt.basis_cross_validate = true; % Maybe always utilize cross-validation in basis tuning : %Turn off cross validation in basis identification. Can save significant computation.
    basis_opt.expand_coeff = 1.1; % Should always use +1 order basis expansion % Order expansion is done multiplicatively (with ceiling function) Lower numbers insure+1 effectiveness.
   % solver_opt.n_folds = 24; % Number of cross-validation folds. Best just as solver_opt.n_folds_pw
    basis_opt.order_bounds = false; % Alternative expansion methods should render this unnecessary.
    
    % Evaluation Related Options
    eval_opt.grad = false; % Requires compatible RHS evaluation. Usually not used.    
    eval_opt.grad_dim = 0; % Number of dimensions with gradients. Zero in this case.
    
    eval_opt.max_dim = dim; % Maximum dimension of problem.
    eval_opt.p_type = repmat(eval_type,1,eval_opt.max_dim);  % Type of polynomial 'l' for legendre, 'h' for hermite. Also defines sample rvs.
    eval_opt.qoi_handle = qoi_handle;

    % Sample Related Options
    sample_opt.sr = 0; % Minimal Sampling Rate
    sample_opt.initial_size = 10; % Useful starting point for many problems.
    sample_opt.n_workers = pool_data.NumWorkers;
    sample_opt.min_sample_percent = 0.1; % Minimum fraction of current sample size to generate at each iteration
    sample_opt.max_sample_percent = 1; % Maximum fraction of current sample size to generate at each iteration
    sample_opt.burn_in = 5000; % Number of burn-in samples
    sample_opt.log_col_rate = -32; % Logarithm of collision rate (MCMC samples to be discarded)
    sample_opt.w_handle = @l2_w; % Minimize l2-coherence
    sample_opt.w_correction_handle = @l2_correction_w; % Minimize l2-coherence with correction sampling
    sample_opt.prop_handle = @orth_prop; % input proposals are from orthogonality measure
    sample_opt.max_samps = 3000; % Maximum number of aggregate samples

    % Solver Related Options
    solver_opt.cross_validate = true; % Enhances stability of estimates and accuracy of error estimate
    solver_opt.n_workers = pool_data.NumWorkers; % Number of parallel workers from pool to use
    solver_opt.n_folds_pw = 6; % Number of CV folds evaluated per worker
    solver_opt.cv_mesh_size = 20; % Number of error tolerances in each cross-validation search.
    solver_opt.solver_handle = @bpdn; % Solver used for regression
    solver_opt.log_tol_exp = 0.1; % Amount for tolerance adjustment self-adjusts well
    solver_opt.comp_samp_perc = 0.8; % Proportion of samples to use for computation, one minus this is proportion for validation
   
    basis_opt.max_elems = 3000; % Maximum number of elements in expansion, generally scales with sample_opt.max_samps
    basis_opt.type_handle = @basis_anisotropic_total_order; % Initial basis is 'to' total order
    basis_opt.dim = dim; % Number of initial dimensions
    basis_opt.ord = ones(1,basis_opt.dim); % Initial total order of approximation
    basis_opt.pc_flag = false; % Whether or not to identify preconditioner, mostly useful for gradient methods
    basis_opt.validation_strikes = 6; % Number of strikes for basis validation
    basis_opt.dim_add = 2; % Number of dimensions to add order 1 approximation to.
    basis_opt.validation_iters = 3; % Number of iterations for basis adaptation.
    basis_opt.basis_identify_handle = @basis_validation_anisotropic_total_order; % Validated basis are anisotropic total order
end
