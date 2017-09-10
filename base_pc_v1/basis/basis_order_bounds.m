% computes order bounds for high-dimensional problems (upper bound only)
% basis_opt = basis_order_bounds(basis_opt, eval_opt)
% -----
% Input
% -----
% basis_opt = options associated with basis
% eval_opt = options related to evaluation
% ------
% Output
% ------
% basis_opt = options associated with basis
function basis_opt = basis_order_bounds(basis_opt, eval_opt)
    [orders,indices] = unique(basis_opt.ord, 'last'); % Here using global minimum order 1 all dimensions
    if isempty(indices) % Sometimes nothing in indices
%        basis_opt.lb = [ones(1,basis_opt.dim_add) zeros(1,eval_opt.max_dim-basis_opt.dim_add)]; % Here using global minimum order 1 all dimensions, otherwise must be computed differently
        basis_opt.ub = basis_opt.lb;
        return
    end
    ub_indices = min(indices+basis_opt.dim_add,eval_opt.max_dim);
    basis_opt.ub = zeros(1,eval_opt.max_dim);
    for kkk = 1:size(orders,2)
        basis_opt.ub(1:ub_indices(kkk)) = max(basis_opt.ub(1:ub_indices(kkk)),orders(kkk));
    end
    basis_opt.ub(1:basis_opt.dim_add) = basis_opt.ub(1:basis_opt.dim_add)+1;
%    o = horzcat(basis_opt.ord,zeros(1,eval_opt.max_dim-basis_opt.dim));
%    basis_opt.lb = max([ones(1,basis_opt.dim_add) zeros(1,eval_opt.max_dim-basis_opt.dim_add)],o-1); % Here using global minimum order 1 all dimensions, otherwise must be computed differently
end
