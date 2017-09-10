% contracts basis set based on solution coefficients
% -----
% basis = basis_coord_max_identify(basis)
% -----
% Input
% -----
% basis = basis object
% ------
% Output
% ------
% basis = basis with corrected .max field
function basis = basis_max_identify(basis) 
basis.max = zeros(1,basis.n_dim); % Must recompute maximum orders and dimension
basis.max_ord = 0;
basis.n_dim = 0;
for k = 1:basis.n_elems
    this_ord = basis.ord{k};
    this_dim = basis.dim{k};
    x = size(this_dim,2);
    s = sum(this_ord);
    if sum(this_ord) > basis.max_ord
        basis.max_ord = s;
    end
    for kk = 1:x
        if(this_ord(kk) > basis.max(this_dim(kk)))
            basis.max(this_dim(kk)) = this_ord(kk);
        end
        if max(this_dim) > basis.n_dim
            basis.n_dim = max(this_dim);            
        end
    end
end
basis.max = basis.max(1:basis.n_dim);
basis.max_ord = max(basis.max);
end
