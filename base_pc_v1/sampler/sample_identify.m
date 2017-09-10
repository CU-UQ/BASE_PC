% updates sample object
% -----
% sample = sample_identify(basis, sample, sample_opt, eval_opt)
% -----
% Input
% -----
% basis = basis object
% sample = sample object
% sample_opt = options for sampling inputs
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% sample = sample object
function sample = sample_identify(basis, sample, sample_opt, eval_opt)
    sample_opt.old_basis = basis; % Old basis needed to correct samples
    sample = sample_adjust(basis, sample, sample_opt, eval_opt); % Adjusts old samples to new basis, applies correction weights
    sample = sample_expand(basis, sample, sample_opt, eval_opt); % Expands samples for new basis sampling, using correction sampling
    sample = folds_identify(sample, solver_opt); % Identification of folds
end
