% applies correction sampling to sample_old
% -----
% sample_new = sample_expand(basis, sample_old, sample_opt, eval_opt)
% -----
% Input
% -----
% basis = basis object
% sample_old = sample object
% sample_opt = options for sampling inputs
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% sample_new = sample object with correction sampling
function sample_new = sample_expand(basis, sample_old, sample_opt, eval_opt) % Expands samples for new basis   
    % generate additional samples for new basis
    sample_new = sample_init(basis, sample_opt, eval_opt); % Need normalizing constant for new sample to adjust several paraemeters
    full_samps = ceil(basis.n_elems*sample_opt.min_sample_percent); % Number of samples after adjustment
    trans_samps = max(ceil(sample_old.n_samps*sample_opt.min_sample_percent),full_samps - sample_old.n_samps); % Number of transition samples is at least a fraction of the number of available samples
    max_samps =  min(sample_opt.max_samps,ceil(sample_old.n_samps*sample_opt.max_sample_percent)); % Maximum number of samples based on number of available samples
    full_samps = sample_old.n_samps + trans_samps; % This is the actual sample size that is kept
    actual_alpha = sample_old.n_samps/full_samps;  % This is the actual alpha, the fraction of old samples
    sample_opt.ac = actual_alpha*sample_old.nc/sample_new.nc; % Adjustment constant (alpha) for adjusting sample
    adj_con = sample_opt.ac; % Will be adjusted as needed
    % Test if further corrections needed
    sample_trans.rej_con = 0; % Initialize value for rej_con (0 implies chain initialization)
    sample_trans.rej_rate = 0; % Initialize value for rej_rate
    sample_trans.n_tries = 0; % Initialize number of tries. Is used to compute rejection rate.
    flag = 1;
    while flag
        [sample_trans, flag] = mcmc_correction_sampler(sample_trans, basis, trans_samps, sample_opt, eval_opt); % If goes through, then all samples are adjusted
        if flag % Sampling not enough to adjust
            sample_opt.ac = flag;
            sample_trans.rej_con = 0; % Reset value
            sample_trans.rej_rate = 0; % Reset value
            sample_trans.n_tries = 0; % Reset number of tries.
            % Adjusts sample size towards maximum or what ac suggests
            if trans_samps < max_samps % No need to check if trans_samps >= max_samps
                prop_rate = (trans_samps*adj_con/sample_opt.ac)/sample_old.n_samps;
                fprintf('Increasing Sample Size. Proposed Sample Rate: %f\n', prop_rate)
                test_val = ceil(trans_samps*adj_con/sample_opt.ac); % Would be our new number of samples
                if test_val >= max_samps;
                    trans_samps = max_samps; % We can only sample up to max_samps
                    full_samps = sample_old.n_samps + trans_samps; % full_samps is adjusted
                    actual_alpha = sample_old.n_samps/full_samps; % actual_alpha is adjusted
                    adj_con = (actual_alpha*sample_old.nc/sample_new.nc); % New adjusting constant has hit its maximum
                    fprintf('Maximum Sample Size Exceeded for Perfect Correction.\n')
                else
                    trans_samps = test_val; % Adjust up to test_val samps
                    full_samps = sample_old.n_samps + trans_samps; % full_samps is adjusted
                    actual_alpha = sample_old.n_samps/full_samps; % actual_alpha is adjusted
                    adj_con = (actual_alpha*sample_old.nc/sample_new.nc); % New adjusting constant for sample
                    sample_opt.ac = adj_con; % because ac is slightly bigger than needed for integerizing.
                end
            end
        end
    end
    adj_con = sample_opt.ac/adj_con; % For weight correction, if ac smaller, than some effective sampling lost
    
    % Concatenate samples 
    sample_new.n_samps = full_samps;
    sample_new.lhs = vertcat(sample_old.lhs, sample_trans.lhs); % Important that sample_old.lhs has been adjusted to new basis
    sample_new.rv = vertcat(sample_old.rv, sample_trans.rv);
    sample_new.p = vertcat(sample_old.p,sample_trans.p); % These numbers will not change
    sample_new.p_rat = vertcat(sample_old.p_rat,sample_trans.p_rat); % These numbers are kept to numbers at the time of sampling
    sample_new.rej_con_history = vertcat(sample_old.rej_con_history,sample_trans.rej_con_history);  % These numbers are kept to numbers at the time of sampling
    sample_new.w = vertcat(adj_con*sample_old.w*sample_old.nc, sample_trans.w)/sample_new.nc; % Important that sample_old.w has been adjusted to new basis sample_old.w already normalized by normalizing constant
    sample_new.wlhs = apply_weights(sample_new.w,sample_new.lhs); % apply weights
    sample_new.rhs = vertcat(sample_old.rhs,eval_opt.qoi_handle(sample_trans.rv,eval_opt)); % Old rhs plus new evaluations
    sample_new.wrhs = apply_weights(sample_new.w,sample_new.rhs); % apply weights
end
