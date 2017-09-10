% Opening pool
delete(gcp('nocreate'))
pool_data = parpool; % Key parts of code are paralellized

% Sample Related Options
sample_opt.sr = 0; % Minimal Sampling Rate
sample_opt.initial_size = 100; % Occasionally want initial_size to differ from sample_opt.sr
sample_opt.min_sample_percent = 0.1; % Percent new samples to generate
sample_opt.n_workers = pool_data.NumWorkers; % Number of parallel workers from pool to use.
sample_opt.burn_in = 25000; % Number of burn-in samples
sample_opt.log_col_rate = -32; % Logarithm of collision rate (MCMC samples to be discarded)
sample_opt.w_handle = @l2_w; % Minimize l2-coherence
sample_opt.w_correction_handle = @l2_correction_w; % Minimize l2-coherence
sample_opt.prop_handle = @orth_prop; % input proposals are from orthogonality measure
sample_opt.max_sample_percent = 0.3; % Maximum proportion of samples to add on each iteration

% Solver Related Options
solver_opt.cross_validate = true;
solver_opt.n_workers = pool_data.NumWorkers; % Number of parallel workers from pool to use
solver_opt.n_folds_pw = 1; % Number of cross-validation folds per worker.
solver_opt.cv_mesh_size = 5; % Number of points in each cross-validation search.
solver_opt.solver_handle = @bpdn; % Solver used for regression
solver_opt.log_tol_exp = 0.1; % Amount for tolerance adjustment
solver_opt.comp_samp_perc = 0.8; % Percentage of samples to use for computation, 1- this is percentage for validation
solver_opt.basis_cross_validate = true; % Turn off cross validation in basis identification. Can save significant computation, when false, but can also lead to unsatisfying basis adaptation.

% Evaluation Related Options
eval_opt.grad = false; % Requires compatible RHS evaluation
eval_opt.max_dim = 3; % Maximum dimension of problem
eval_opt.grad_dim = 0; % Number of dimensions with gradients. (must correspond to first dimensions)
eval_opt.p_type = repmat('h',1,eval_opt.max_dim);  % type of polynomial 'l' for legendre, 'h' for hermite. Also defines sample rvs.
eval_opt.qoi_handle = @duffing_hermite_eval;
eval_opt.parallel_qoi_eval = true; % Evaluates multiple QoI at one time Useful when QoI Eval is expensive.

% Convergence Options
max_iter = 100; % Maximum sample_adaptation iterations
basis_opt.max_elems = 3000; % Maximum number of elements in expansion
sample_opt.max_samps = 3000; % Maximum number of samples to draw

% ---------- %
% initialize %
% ---------- %

for kk = 1:5
    % Basis Related Options can be changed on iteration (these are just
    % repeated from above
    basis_opt.type_handle = @basis_anisotropic_total_order; % Initial basis is 'to' total order
    basis_opt.dim = 3; % Number of initial dimensions
    basis_opt.ord = 3*ones(1,basis_opt.dim); % Initial total order of approximation
    basis_opt.pc_flag = false; % Whether or not to identify preconditioner
    basis_opt.validation_strikes = 6; % Number of strikes for basis validation
    basis_opt.expand_coeff = 1.1; % Order expansion is done multiplicatively (with ceiling function) Lower numbers insure+1 effectiveness.
    basis_opt.dim_add = 3; % Number of dimensions to add order 1 approximation to.
    basis_opt.validation_iters = 3; % Number of iterations for basis adaptation
    basis_opt.order_bounds = false;
    basis_opt.basis_identify_handle = @basis_validation_anisotropic_total_order;
    
    basis = basis_init(basis_opt, eval_opt); % Initialization of basis    
    sample = sample_init(basis, sample_opt, eval_opt); % Initialization of Sample
    sample = folds_identify(sample, solver_opt); % Identification of folds
    solver_opt.sig = 1; % Initialize error estimate
    [surrogate, solver_opt] = solution_identify(sample,solver_opt); % generates surrogate
    error = []; %surrogate.err; % For saving for diagnostics
    tol = []; %surrogate.sig; % For diagnostics, identified via cross-validation
    n_elems = []; % basis.n_elems; % For diagnostics: Number of elements in basis
    n_evals = []; %sample.n_samps; % For diagnostics: Number of evaluations
    ord = []; % For diagnostics: Holds order information
    dim = [];
    % ord{1} = basis.max; % Adding initialization.
    fprintf('Initialization. Samples: %i  Basis Elements: %i Validated Error: %f Tolerance: %f\n', sample.n_samps, basis.n_elems, surrogate.err, surrogate.sig);
    
% ---- %
% loop %
% ---- %
    k = 0;
    while(k < max_iter && sample.n_samps < sample_opt.max_samps) % Several exit criteria    
        sample_opt.old_basis = basis; % Old basis needed to correct samples
        [basis, basis_opt, solver_opt, sample, ~] = basis_identify(basis, sample, surrogate, basis_opt, eval_opt, solver_opt);
        sample = sample_adjust(basis, sample, sample_opt, eval_opt); % Adjusts old samples to new basis
        sample = sample_expand(basis, sample, sample_opt, eval_opt); % Expands samples for new basis sampling, using correction sampling
        sample = folds_identify(sample, solver_opt); % Identification of folds
        [surrogate, solver_opt] = solution_identify(sample, solver_opt); % Updates surrogate
        error = [error; surrogate.err]; %#ok<AGROW>
        tol = [tol; surrogate.sig]; %#ok<AGROW>
        n_elems = [n_elems; basis.n_elems]; %#ok<AGROW>
        n_evals = [n_evals; sample.n_samps]; %#ok<AGROW>
        k = k+1;
        ord{k} = basis.max; %#ok<SAGROW>
        dim = [dim; basis.n_dim]; %#ok<AGROW>
        fprintf('Iteration: %i. Samples: %i Basis Elements: %i Validated Error: %f Tolerance: %f\n', k, sample.n_samps, basis.n_elems, surrogate.err, surrogate.sig);
    str = strcat('duffing_ex_sa_',num2str(kk));
    save(str)
    end
    
% ------ %
% finish %
% ------ %
    save(str)
end
