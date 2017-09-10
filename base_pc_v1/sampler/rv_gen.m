% Generates random variable from orthogonality distribution (no importance sampling)
% -----
% [rv, p] = rv_gen(opt)
% -----
% Input
% -----
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% rv = random vector
% p = probability of draw
function [rv,p] = rv_gen(eval_opt)
    rv = zeros(1,eval_opt.max_dim);
    p = 1;
    for kk = 1: eval_opt.max_dim % Generate test random variable
        switch eval_opt.p_type(kk)
            case 'H' % hermite/normal with scale/location
                rv(kk) = randn;
                p = p*normpdf(rv(kk));
                rv(kk) = rv(kk)*eval_opt.sig(kk) + eval_opt.mu(kk);
            case 'h' % hermite/normal
                rv(kk) = randn;
                p = p*normpdf(rv(kk));
            case 'l' % Legendre for uniform [0,1]
                rv(kk) = rand;
            case 'L' % Legendre for uniform [-1,1]
                rv(kk) = 2*rand-1;
                p = 1/2*p; % Uniform rv on [-1,1]
            case 'a' % laguerre/gamma
                rv(kk) = gamrnd(eval_opt.alpha(kk),1);
                p = p*gampdf(rv(kk),eval_opt.alpha(kk),1);
            case 'j' % jacobi/beta [-1,1]
                rv(kk) = betarnd(eval_opt.beta(kk)+1,eval_opt.alpha(kk)+1); % Not yet stretched to [-1, 1]
                p = p*0.5*betapdf(rv(kk),eval_opt.alpha(kk),eval_opt.beta(kk)); % Half probability due to stretch
                rv(kk) = 2*rv(kk)-1; % Jacobi rv we place in [-1,1] here, as opposed to the legendre rv.
        end
    end
end
