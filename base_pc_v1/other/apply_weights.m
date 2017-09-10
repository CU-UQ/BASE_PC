% applies weights to matrix
% -----
% mat = apply_weights(w,mat);
% -----
% Input
% -----
% w = weights
% mat = matrix to apply weights to
% ------
% Output
% ------
% mat = matrix with weights applied
function mat = apply_weights(w,mat)
    n1 = size(w,1);
    n2 = size(mat,1);
    k = n2/n1;
    if k == 1
        mat = bsxfun(@times, w, mat);
        return
    end
    set = 1:n1;
    for j = 1:k
        mat(set,:) = bsxfun(@times, w, mat(set,:));
        set = set+1;
    end
end
