% decomposes variance by coordinate using pce coefficients
% -----
% v = variance_by_coord(basis,c)
% -----
% Input
% -----
% basis = basis object
% c = solution coefficients
% ------
% Output
% ------
% v = fraction of variance in each coordiante
function v = variance_by_coord(basis,c)
    v = zeros(1,basis.n_dim);
    tot_v = 0;
    for j = 1:basis.n_elems
        this_v = c(j)^2;
        this_dim = basis.dim{j};
        this_ord = basis.ord{j};
        n_dims = size(this_dim,2);
        if(~isempty(basis.dim{j}))
            tot_ord = sum(this_ord);
            tot_v = tot_v + this_v;
            for i = 1:n_dims
                v(this_dim(i)) = v(this_dim(i)) + this_v*this_ord(i)/tot_ord;
            end
        end
    end
    v = v/tot_v;
end
