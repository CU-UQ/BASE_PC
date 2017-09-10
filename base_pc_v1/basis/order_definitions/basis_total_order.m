% returns total order basis object
% -----------
% basis = basis_total_order(basis_opt)
% -----------
% Input
% -----
% basis_opt = options to describe basis
% ------
% Output
% ------
% basis = basis object
function basis = basis_total_order(basis_opt) % Can be made more efficient through vectorization
    basis.n_dim = basis_opt.dim;
    basis.max = basis_opt.ord.*ones(1,basis_opt.dim);
    basis.max_ord = basis_opt.ord;
    switch basis_opt.ord % A few easy cases
        case 0
            basis.ord={[]};
            basis.dim={[]};
            basis.n_elems = 1;
            return
        case 1            
            basis.dim = [{[]}; num2cell((1:basis_opt.dim)')];
            basis.ord = [{[]}; num2cell(ones(basis_opt.dim,1))];
            basis.n_elems = basis_opt.dim+1;
            return
    end
    if basis_opt.dim == 1 % Another easy cases
        basis.dim = [{[]}; num2cell(ones(basis_opt.ord,1))];
        basis.ord = [{[]}; num2cell((1:basis_opt.ord)')];
        basis.n_elems = basis_opt.ord+1;
        return
    end
    d_minus = basis_opt.dim-1;
    c_test = 1:d_minus;
    % General case
    basis.n_elems = round(exp(gammaln(basis_opt.ord+basis_opt.dim+1)-gammaln(basis_opt.ord+1)-gammaln(basis_opt.dim+1)));
    basis.ord = cell(basis.n_elems,1);
    basis.dim = cell(basis.n_elems,1);
    basis_vec = zeros(1,basis_opt.dim);
    n = 1;
    t = zeros(1,d_minus);
    for c_ord = 1:basis_opt.ord
        end_gen = 0;
        first_this_order = 0;
        c_test = c_test + 1;
        while (end_gen == 0)
            n = n+1;
            if (first_this_order == 0)
                t = 1:d_minus;
                first_this_order = 1;
            else
                l = find(t < c_test ,1,'last'); % Find index to increment
                t(l:d_minus) = (t(l)+1):(t(l)+basis_opt.dim-l); % Increment t(l) and  all t(k) for k > l to be t(l) + k-l. 
            end
            basis_vec(1) = t(1)-1;  
            basis_vec(2:d_minus) = t(2:d_minus) - t(1:(d_minus-1)) - 1; % Subsequent orders are t(k) - t(k-1) - 1 Do not need to check higher dimensions until c_ord increased
            basis_vec(basis_opt.dim) = basis_opt.dim+c_ord-t(d_minus)-1;
            if (t(1) == c_ord+1)
                end_gen = end_gen + 1;
            end
            act_opt = find(basis_vec~=0); % Puts into our basis definition
            basis.dim{n} = act_opt; % Active dimensions
            basis.ord{n} = basis_vec(act_opt); % Active orders
        end
    end
end
