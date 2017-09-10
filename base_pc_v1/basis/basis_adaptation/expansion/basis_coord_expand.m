% expands basis set using given parameters
% -----
% [basis_exp, basis_opt] = basis_coord_expand(basis_c, basis_opt, eval_opt)
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
function [basis_exp, basis_opt] = basis_coord_expand(basis_c, basis_opt, eval_opt, ord_indices)
    basis_opt.ord = max(1,basis_c.max+ord_indices); % Increase order by 1 in new coordinate dimension
    basis_opt.dim = min(basis_c.n_dim+basis_opt.dim_add,eval_opt.max_dim); % In case of adaptive dimensionality
    basis_opt.ord = horzcat(basis_opt.ord, 1*ones(1,basis_opt.dim-size(basis_opt.ord,2))); % Increase dim if needed
    if isfield(basis_c,'hyp') % Not yet implemented well
        basis_opt.hyp = basis_c.hyp; % Hyperbolicity not adjusted
    end
    basis_exp = basis_init(basis_opt, eval_opt);   
end
