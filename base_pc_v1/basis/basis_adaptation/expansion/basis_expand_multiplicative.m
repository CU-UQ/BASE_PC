% expands basis set using given parameters
% -----
% [basis_exp, basis_opt] = basis_expand(basis_c, basis_opt, eval_opt)
% -----
% Input
% -----
% basis_c = contracted basis
% basis_opt = basis options
% eval_opt = options related to evaluation
% ------
% Output
% ------
% basis_exp = expanded basis object
% basis_opt = updated basis options
function [basis_exp, basis_opt] = basis_expand_multiplicative(basis_c, basis_opt, eval_opt)
    basis_opt.dim = min(basis_c.n_dim+basis_opt.dim_add,eval_opt.max_dim); % In case of adaptive dimensionality
    ord = ceil(max(1,basis_c.max*basis_opt.expand_coeff)); % Increase order by at least 1 in each already on dimension (depending on present order), only by 1 in new dimensions 
    % This number can be tuned
    basis_opt.ord = horzcat(ord, ones(1,basis_opt.dim-size(ord,2))); % Increase dim if needed
    if isfield(basis_c,'hyp')
        basis_opt.hyp = basis_c.hyp; % Hyperbolicity not adjusted
    end
    if basis_opt.order_bounds
        basis_opt.ord = horzcat(basis_opt.ord,zeros(1,size(basis_opt.ub,2)-size(basis_opt.ord,2))); % incase basis_opt.ord has dropped some dimensions
        basis_opt.ord = min(basis_opt.ord,basis_opt.ub);
%        basis_opt.ord = max(basis_opt.ord,basis_opt.lb);
        basis_opt.dim = find(basis_opt.ord>0,1,'last'); % Can change when using bounds.
        basis_opt.ord = basis_opt.ord(1:basis_opt.dim);
    end
    basis_exp = basis_init(basis_opt, eval_opt);
end
