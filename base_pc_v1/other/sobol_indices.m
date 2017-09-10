% Computes Sobol Indices and estimate variance
% -----
% [sobol, dim_index var] = sobol_indices(c, basis)
% -----
% Input
% -----
% c = coefficients associated with basis
% basis = associated basis
% ------
% Output
% ------
% sobol = Sobol indices sorted in decreasing order
% dim_index = dimensions sorted in decreasing order of sobol index
% var = total estimated variance

function [sobol, dim_index, var] = sobol_indices(c, basis)

var = 0; % Will be computed as well

sobol = zeros(basis.n_dim,1);
for i = 1:basis.n_elems
    if ~isempty(basis.ord{i})
        var = var + c(i).^2;
        for k = 1:size(basis.dim{i},2)
            sobol(basis.dim{i}(k)) = sobol(basis.dim{i}(k)) + c(i).^2;
        end
    end
end

sobol = sobol/var;

[sobol,dim_index] = sort(sobol,'descend');

end