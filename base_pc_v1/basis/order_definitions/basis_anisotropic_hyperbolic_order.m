% returns anisotropic, hyperbolic, total order basis object
% -----------
% basis = basis_anisotropic_hyperbolic_order(basis_opt)
% -----------
% Input
% -----
% basis_opt = options to describe basis
% ------
% Output
% ------
% basis = basis object
function basis = basis_anisotropic_hyperbolic_order(basis_opt)
    % total order is somewhat analogous to the sum of orders
    max_ord = max(basis_opt.ord);
    basis.max_ord = max(basis_opt.ord);
    basis.max = basis_opt.ord;
    basis.hyp = basis_opt.hyp;
    basis.n_dim = basis_opt.dim;
    d_minus = basis_opt.dim-1;
    n = 1;
    t = zeros(1,d_minus);
    basis.ord={[]};
    basis.dim={[]};
    % The basis functions are generated and tested

    for c_ord = 1:max_ord % Loop over potential basis functions
        first_this = 0;
        basis_array_prop = zeros(1,basis_opt.dim);
        while true % break loop if condition follows
            if (first_this == 0)
                for k=1:d_minus
                    t(k) = k;
                end
                first_this = 1;
            else
                if (t(d_minus) < d_minus + c_ord)
                    t(d_minus) = t(d_minus) + 1;
                else
                    l=d_minus;
                    while (t(l) == l+c_ord)
                        l=l-1;
                    end
                    t(l) = t(l)+1;
                    for k =l+1:d_minus
                        t(k) = t(l)+k-l;
                    end
                end
            end
            basis_array_prop(1) = t(1)-1;
            for k= 2:d_minus
                basis_array_prop(k) = t(k)-t(k-1)-1;
            end
            basis_array_prop(basis_opt.dim) = basis_opt.dim+c_ord-t(d_minus)-1; % Final entry treated a bit differently
            if(sum((basis_array_prop./basis_opt.ord).^basis_opt.hyp) <= 1)
                n = n+1;
                act_dim = find(basis_array_prop~=0);
                act_ord = basis_array_prop(act_dim);
                basis.dim{n} = act_dim;
                basis.ord{n} = act_ord;
                if sum(act_ord) > basis.max_ord;
                    basis.max_ord = sum(act_ord);
                end
            end
            if (t(1) == c_ord+1)
                break % break loop
            end
        end
    end
    basis.dim = basis.dim';
    basis.ord = basis.ord';
    basis.n_elems = n;
end
