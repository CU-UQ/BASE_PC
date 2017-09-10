delete(gcp('nocreate'))
pool_data = parpool; % Key parts of code are paralellized
dim = 2;
eval_type = 'l';
qoi_handle = @exp_sine_decay_eval;
[basis_opt, eval_opt, sample_opt, solver_opt] = default_parameters(pool_data,dim,eval_type,qoi_handle);
eval_opt.exp_const = 1; % Used with @exp_decay_eval
eval_opt.exp_decay_init = 2; % Used with @exp_decay_eval
eval_opt.exp_decay_exp = 1; % used with exp_decay_eval

basis = basis_init(basis_opt, eval_opt); % Identifies basis.
fprintf('Basis Initialized\n');
sample = sample_init(basis, sample_opt, eval_opt);
fprintf('Sample Initialized\n');
sol = solution_identify(sample,solver_opt); % Identifies solution.
fprintf('Cross Validated Solution Computed\n');