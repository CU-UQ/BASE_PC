delete(gcp('nocreate'))
pool_data = parpool; % Key parts of code are paralellized
% Evaluation Related Options
eval_opt.qoi_handle = @exp_sine_decay_eval; % Handle for evaluating QoI
    eval_opt.exp_const = 1; % Used with @exp_decay_eval
    eval_opt.exp_decay_init = 2; % Used with @exp_decay_eval
    eval_opt.exp_decay_exp = 1; % used with exp_decay_eval
eval_opt.grad = false; % Whether or not gradient terms are evaluated. Requires corresponding QoI handle.
eval_opt.max_dim = 2; % Maximum dimension of problem
eval_opt.grad_dim = 0; % Number of dimensions with gradients. (must correspond to first dimensions used)
eval_opt.p_type = repmat('L',1,eval_opt.max_dim);  % Type of polynomial 'L' here means legendre in each dimension. Also used to determine random variables used.
eval_opt.parallel_qoi_eval = false; % Evaluates multiple QoI at one time Useful when QoI Eval is expensive.

% Basis Related Options
basis_opt.type_handle = @basis_total_order; % Initial basis is total order
basis_opt.dim = eval_opt.max_dim; % Number of initial dimensions. Here, we use all dimensions.
basis_opt.ord = 2; % Initial total order of approximation. Here two.
basis_opt.pc_flag = false; % Whether or not to identify preconditioner. This option is useful for research if you define such a preconditioner.

% Options for generating a candidate sample set
sample_opt.sr = 0; % Can maintain a ratio of samples to basis functions if wanted, rarely useful for l1
sample_opt.initial_size = 15000; % Initial size of sample.
sample_opt.n_workers = pool_data.NumWorkers; % Number of parallel workers from pool to use for sampling.
% These determine what distribution samples are drawn from
sample_opt.w_handle = @one_w; % Leads to sampling from orthogonality distribution.
sample_opt.prop_handle = @orth_prop; % Proposals from orthogonality distribution, so sampling is exact.
% These are for quality of MCMC sampler distribution matching target
sample_opt.burn_in = 0; % Number of burn-in samples, unnecessary for orthogonal distribution sampling
sample_opt.log_col_rate = -1; % Logarithm of collision rate, unnecessary for orthogonal distribution sampling

% Solver Related Options
solver_opt.cross_validate = true;
solver_opt.n_workers = pool_data.NumWorkers; % Number of parallel workers from pool to use
solver_opt.tol = 0; % Stop if cross-validated fraction-of-variance-unexplained is less than this.
solver_opt.n_folds_pw = 2; % Number of cross-validation folds per worker. More folds implies more accurate cross-validation.
solver_opt.cv_mesh_size = 20; % Number of points in each cross-validation search.
solver_opt.solver_handle = @bpdn; % Solver used for regression
solver_opt.log_tol_exp = 0.3; % Parameter for determining how largest cross-validation candidate values are generated.
solver_opt.comp_samp_perc = 0.8; % Percentage of samples per computation. Remaining samples are used for validation.
solver_opt.err = 1; % Without any computed solution, initial error is 1 i.e. 100%

basis = basis_init(basis_opt, eval_opt); % Identifies basis.
fprintf('Basis Initialized\n');
sample = orth_sampler(sample_opt, eval_opt); % Identifies sample.
sample.wlhs = basis_eval(basis,sample.rv,eval_opt);
sample.wrhs = sample.rhs;
fprintf('Sample Initialized\n');
sol = solution_identify(sample,solver_opt); % Identifies solution.
fprintf('Cross Validated Solution Computed\n');
