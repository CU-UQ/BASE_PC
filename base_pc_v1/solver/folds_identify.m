% identifies new validation folds when sample size is increased
% -----
% sol = folds_identify(sample,solver_opt)
% -----
% Input
% -----
% sample = sample object
% solver_opt = options for identifying solution
% ------
% Output
% ------
% sample = sample object
function sample = folds_identify(sample, solver_opt)
    sample.n_folds = solver_opt.n_workers*solver_opt.n_folds_pw;    
    sample.folds = zeros(sample.n_folds,sample.n_samps);
    for k = 1:sample.n_folds % Build fold sets
        sample.folds(k,:) = randperm(sample.n_samps);
    end
end
