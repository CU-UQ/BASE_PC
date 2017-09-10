% proposal distribution based on asymptotically motivated distributions depending on inputs
% -----
% [rv,p,p_rat] = asym_prop(sample_opt, eval_opt)
% -----
% Input
% -----
% sample_opt = options for sampling inputs
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% rv = sampled inputs
% p = probability of sample in orthogonality distribution. Used for diagnostics
% p_rat = ratio of proposal distribution over orthogonality distribution. Used for MCMC sampler
function [rv,p,p_rat] = asym_prop(sample_opt, eval_opt)
    rv = zeros(1,eval_opt.max_dim);
    p = 1;
    p_rat = 1;
    for k = 1: eval_opt.max_dim %% Generate test random variable
        switch eval_opt.p_type(k)
            case 'h' % hermite/normal Asymp approx is uniform (Note this is individual dimension `asymptotic', and not even true asymptotic)
                r = sqrt(2)*sqrt(2*sample_opt.order(k)+2); % radius, provides diameter of uniform sampling
                rv(k) = r*(2*rand-1);
                temp = normpdf(rv(k));
                p = p*temp;
                p_rat = p_rat*(2*r)/temp;
            case 'l' % legendre/uniform [0,1] Asymp is shifted Cheby
                rv(k) = 0.5*cos(pi*rand)+0.5; % Rescale so rv is in [0,1]
                %p = p; % pdf = 1
                p_rat = p_rat*2/(pi*(1-(2*rv(k)-1)^2)^(0.5)); % Appropriate ratio for rv in [0,1]
            case 'L' % legendre/uniform [-1,1] Asymp is cheby
                rv(k) = cos(pi*rand); % rv in [-1,1]
                p = 1/2*p;
                p_rat = p_rat*2/(pi*(1-rv(k)^2)^(0.5)); % Appropriate ratio for rv in [-1,1]
            case 'a' % laguerre/gamma Asymp is different gamma
                rv(k) = gamrnd(eval_opt.alpha(k)+2,1);
                temp = gampdf(rv(k),eval_opt.alpha(k)+1,1);
                p = p*temp;
                p_rat = p_rat*gampdf(rv(k),sample_opt.alpha(k)+2,1)/temp;
            case 'j' % jacobi/beta [-1,1] Asymp is different beta
                rv(k) = betarnd(eval_opt.beta(k)+0.5,eval_opt.alpha(k)+0.5); % Not yet stretched to [-1, 1]
                temp = 0.5*betapdf(rv(k),eval_opt.beta(k)+0.5,eval_opt.alpha(k)+0.5);
                p = p*temp; % Half probability due to stretch
                p_rat = p_rat*0.5*betapdf(rv(k),eval_opt.beta(k)+1,eval_opt.alpha(k)+1)/temp;
                rv(k) = 2*rv(k)-1; % Jacobi rv in [-1,1]
        end
    end
end
