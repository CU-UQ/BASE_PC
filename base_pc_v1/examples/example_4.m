delete(gcp('nocreate'))
pool_data = parpool; % Key parts of code are paralellized
dim = 2;
eval_type = 'l';
qoi_handle = @franke_eval;
[basis_opt, eval_opt, sample_opt, solver_opt] = default_parameters(pool_data,dim,eval_type,qoi_handle);
sample_opt.w_handle = @l2_w; % Leads to sampling from orthogonality distribution.
sample_opt.prop_handle = @orth_prop; % Proposals from orthogonality distribution, so sampling is exact.
sample_opt.w_correction_handle = @l2_correction_w; % Minimize l2-coherence with correction sampling
  

% eval_opt.exp_const = 1; % Used with @exp_decay_eval
% eval_opt.exp_decay_init = 2; % Used with @exp_decay_eval
% eval_opt.exp_decay_exp = 1; % used with exp_decay_eval

basis = basis_init(basis_opt, eval_opt); % Identifies basis.
fprintf('Basis Initialized\n');
sample = sample_init(basis, sample_opt, eval_opt);
fprintf('Sample Initialized\n');
sol = solution_identify(sample,solver_opt); % Identifies solution.
fprintf('Cross Validated Solution Computed\n');

basis_opt.basis_identify_handle = @basis_validation_anisotropic_total_order;
k = 1;
sample_opt.old_basis = basis;
max_iter = 50;
while(k < max_iter && sample.n_samps < sample_opt.max_samps) % Several exit criteria
    % Begins with basis adaptation
    sample = sample_expand(basis, sample, sample_opt, eval_opt); % Expands samples for new basis sampling, using correction sampling
    fprintf('Samples Expanded.\n')
    sample = folds_identify(sample, solver_opt); % Identification of folds
    fprintf('Folds Identified.\n');
    [surrogate, solver_opt] = solution_identify(sample, solver_opt); % Updates surrogate
    fprintf('Initial Solution Identified.\n')
    sample_opt.old_basis = basis; % Old basis needed to correct samples
    [basis, basis_opt, solver_opt, sample, surrogate] = basis_identify(basis, sample, surrogate, basis_opt, eval_opt, solver_opt); % All adaptivity cuz of basis adaptation
    fprintf('Adjusted Basis Identified.\n')
    k = k+1;
end