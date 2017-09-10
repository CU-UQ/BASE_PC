% returns anisotropic total order basis object
% -----------
% basis = basis_anisotropic_total_order(basis_opt)
% -----------
% Input
% -----
% basis_opt = options to describe basis
% ------
% Output
% ------
% basis = basis object
function basis = basis_anisotropic_total_order(basis_opt)
    % total order is somewhat analogous to the sum of orders
    max_ord = max(basis_opt.ord);
    basis.max = basis_opt.ord;
    basis.max_ord = max(basis_opt.ord);
    basis.n_dim = basis_opt.dim;
    if basis.n_dim == 1 % We only allow dimensions to be defined from first. This one is done.
        basis.ord={[]};
        basis.dim={[]};
        for k = 1:basis_opt.ord(1)        
            basis.dim{k+1} = 1;
            basis.ord{k+1} = k;
        end
        basis.dim = basis.dim';
        basis.ord = basis.ord';
        basis.n_elems = basis_opt.ord(1)+1;
        return
    end
    d_minus = basis_opt.dim-1; % we use this dim-1 repeatedly
    n = 1;
    t = zeros(1,d_minus); % loops over potential basis functions
    basis.ord={[]};
    basis.dim={[]};
    c_test = 1:d_minus;

    [m_ord, index] = sort(basis_opt.ord, 'ascend'); % We sort this ascending to determine break point.
    [~, index] = sort(index); % This is the permutation indices we require
    % The basis functions are generated and tested

    for c_ord = 1:max_ord % Loop over potential basis functions
        c_test = c_test + 1; % Increment c_test
        first_this = 0;
        basis_array_prop = zeros(1,basis_opt.dim);
        while true % break loop if condition follows
            if (first_this == 0)
                t = 1:d_minus; % Corresponds to many zero orders
                first_this = 1; % After initialization no longer needed
            else
                l = find(t < c_test ,1,'last'); % Find index to increment
                t(l:d_minus) = (t(l)+1):(t(l)+basis_opt.dim-l); % Increment t(l) and  all t(k) for k > l to be t(l) + k-l. 
                    % can perhaps accelerate this more by redefining d_minus appropriately. Many t(k) are redefined to same value
            end
            basis_array_prop(1) = t(1)-1; % First order is t(1)-1
            basis_array_prop(2:d_minus) = t(2:d_minus) - t(1:(d_minus-1)) - 1; % Subsequent orders are t(k) - t(k-1) - 1 Do not need to check higher dimensions until c_ord increased
            basis_array_prop(basis_opt.dim) = basis_opt.dim+c_ord-t(d_minus)-1; % Final entry fills in gap to c_ord
            if(sum(basis_array_prop./m_ord) <= 1) % If this condition, then valid basis function
                n = n+1;
                basis_array_prop = basis_array_prop(index);
                act_dim = find(basis_array_prop~=0);
                act_ord = basis_array_prop(basis_array_prop~=0); % Order for active dimensions
                basis.dim{n} = act_dim; % How basis dimension is stored
                basis.ord{n} = act_ord; % How basis order is stored
            else % If condition not held then can safely ignore a number of other vectors due to ascended sorting
                if(l == 1) % If on first dimension, no others fit
                    break
                else % Can safely remove any other incrementation in this dimension
                    t(l:d_minus) = c_test(l:d_minus);
                end
            end
            if (t(1) >= c_ord+1) % basis search exhausted for this c_ord
                break % break loop
            end
        end
    end
    basis.dim = basis.dim';
    basis.ord = basis.ord';
    basis.n_elems = n;
end
