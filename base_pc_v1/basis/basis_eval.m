% evaluates basis functions at input
% -----
% lhs = basis_eval(basis, input)
% -----
% Input
% -----
% basis = basis object
% input = input evaluated according to basis
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% lhs = matrix of basis functions evaluated at input
function lhs = basis_eval(basis,input,eval_opt)
    n_samps = size(input,1);
    n1 = min(eval_opt.grad_dim,basis.n_dim)+1;
    n_rows = n1*n_samps;
    lhs = ones(n_rows,basis.n_elems);
    for k=1:basis.n_dim % iterates over dimension
        % May be more prudent iteration over basis functions
        if floor(basis.max(k)) > 0 % nothing to do if no effective order
            if eval_opt.grad && k <= eval_opt.grad_dim
                f_set = 1:n1:n_rows;
                fp_set = (k+1):n1:n_rows;
                switch eval_opt.p_type(k)
                    case 'H' % Hermite/ normal with location/scale parameter
                        ord = floor(basis.max(k));
                        [f, fp] = hermite_eval_gradient(ord,(input(:,k)-eval_opt.mu(k))/eval_opt.sig(k));
                        nc = 1;
                        for kk = 1:ord
                            nc = nc*eval_opt.sig(k);
                            f(:,kk) = f(:,kk)/nc;
                            fp(:,kk) = f(:,kk)/nc;
                        end % Repositions/Rescales hermite polynomials
                    case 'h' % Hermite/ normal                    
                        [f,fp] = hermite_eval_gradient(floor(basis.max(k)),input(:,k));
                    case 'l' % Legendre/ Uniform [0,1]
                        [f,fp] = legendre_eval_gradient(floor(basis.max(k)),2*input(:,k)-1); % Adjusts for differences between poly orthogonality and rv gen
                    case 'L' % Legendre/ Uniform [-1,1]
                        [f,fp] = legendre_eval_gradient(floor(basis.max(k)),input(:,k)); % Adjusts for differences between poly orthogonality and rv gen
                    case 'a' % Laguerre / Gamma
                        [f,fp] = laguerre_eval_gradient(floor(basis.max(k)),input(:,k),eval_opt.alpha(k)-1); % Adjusts for differences between poly orthogonality and rv gen
                    case 'j' % Jacobi / Beta [-1,1]
                        [f,fp] = jacobi_eval_gradient(basis.max(k),input(:,k),eval_opt.beta(k)-1,eval_opt.alpha(k)-1); % Adjusts for differences between poly orthogonality and rv gen
                end
            else
                switch eval_opt.p_type(k)
                    case 'H' % Hermite/ normal with location/scale parameter
                        ord = floor(basis.max(k));
                        f = hermite_eval(ord,(input(:,k)-eval_opt.mu(k))/eval_opt.sig(k)); % Repositions/Rescales hermite polynomials
                        nc = 1;
                        for kk = 1:ord
                            nc = nc*eval_opt.sig(k);
                            f(:,kk) = f(:,kk)/nc;
                        end
                    case 'h' % Hermite/ normal
                        f = hermite_eval(floor(basis.max(k)),input(:,k));
                    case 'l' % Legendre/ Uniform [0,1]
                        f = legendre_eval(floor(basis.max(k)),2*input(:,k)-1); % Adjusts for differences between poly/rv
                    case 'L' % Legendre/ Uniform [-1,1]
                        f = legendre_eval(floor(basis.max(k)),input(:,k)); % Adjusts for differences between poly/rv
                    case 'a' % Laguerre / Gamma
                        f = laguerre_eval(floor(basis.max(k)),input(:,k),eval_opt.alpha(k)-1); % Adjusts for differences between poly/rv
                    case 'j' % Jacobi / Beta [-1,1]
                        f = jacobi_eval(floor(basis.max(k)),input(:,k),eval_opt.beta(k)-1,eval_opt.alpha(k)-1); % Adjusts for differences between poly/rv
                end
            end
        end
        for kk = 1:basis.n_elems % iterates over basis functions
            l = find(basis.dim{kk} == k,1); % find appropriate dimension index
            if eval_opt.grad && k <= eval_opt.grad_dim
                if isempty(l)
                    lhs(fp_set,kk) = zeros(n_samps,1); % derivative of constant is zero
                else
                    lhs(f_set,kk) = lhs(f_set,kk).*f(:,basis.ord{kk}(l));
                    lhs(fp_set,kk) = lhs(fp_set,kk).*fp(:,basis.ord{kk}(l));
                end
            else
                if ~isempty(l)
                    lhs(:,kk) = lhs(:,kk).*f(:,basis.ord{kk}(l));
                end
            end
        end
    end
    if basis.pc_flag
        lhs = lhs*basis.pc;
    end
end
