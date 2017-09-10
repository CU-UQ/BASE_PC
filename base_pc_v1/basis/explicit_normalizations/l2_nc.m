% normalization constant for l2-coh-opt sampling
% -----------
% c = l2_nc(basis,o)
% -----------
% Input
% -----
% basis = basis_object
% o = option parameters, not used here
% ------
% Output
% ------
% c = normalizing constant

function c = l2_nc(basis,o)
 c = sqrt(basis.n_elems);
end
