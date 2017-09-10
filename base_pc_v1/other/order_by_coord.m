% percentage of basis order in each coordinate
% -----
% v = order_by_coord(basis)
% -----
% Input
% -----
% basis = basis object
% ------
% Output
% ------
% v = fraction of basis order in each coordinate
function v = order_by_coord(basis)
    v = zeros(1,basis.n_dim);
    tot_v = 0;
    for j = 1:basis.n_elems
        this_dim = basis.dim{j};
        this_ord = basis.ord{j};
        n_dims = size(this_dim,2);
        if(~isempty(basis.dim{j}))
            tot_v = tot_v+1;
            tot_ord = sum(this_ord);
            for i = 1:n_dims
                v(this_dim(i)) = v(this_dim(i)) + this_ord(i)/tot_ord;
            end
        end
    end
    v = v/tot_v;
end
