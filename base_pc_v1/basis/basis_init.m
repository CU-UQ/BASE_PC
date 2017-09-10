% constructs initial basis opbject
% -----
% basis = basis_init(basis_opt, eval_opt)
% -----
% Input
% -----
% basis_opt = options to identify basis
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% basis = basis object
function basis = basis_init(basis_opt, eval_opt) % can be modified for different bases
    % Identifies basis based on ib_handle
    basis = basis_opt.type_handle(basis_opt);    
    % If preconditioning is applied
    if basis_opt.pc_flag
        if isfield(basis_opt, 'pc_handle')
            basis.pc = basis_opt.pc_handle(basis, basis_opt);
        else % The below process is slow and perhaps inaccurate, and should probably be removed.
            n_samps = min(10^4,max(10,basis_opt.sr)*basis.n_elems);
            rv = zeros(n_samps,basis_opt.max_dim);
            basis.pc_flag = false;
            for k = 1:n_samps
                rv(k,:) = rv_gen(basis_opt);
            end
            mat = basis_eval(basis,rv, basis_opt, eval_opt);
            mat = mat'*mat/n_samps;
            for kk = 2:basis_opt.burn_in
                parfor k = 1:n_samps
                    rv(k,:) = rv_gen(eval_opt);
                end
                mat1 = basis_eval(basis,rv, basis_opt, eval_opt);
                mat1 = mat1'*mat1/n_samps;
                mat = mat*(kk-1)/kk + mat1/kk;
                basis.pc = inv(chol(mat));
            end
        end        
        basis.pc_flag = true;
    else
        basis.pc_flag = false;
    end
end
