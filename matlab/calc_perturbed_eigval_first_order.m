function [E1 ] = calc_perturbed_eigval_first_order(Qp,vals,vecs,eig_ind)
%[vecs vals] = eig(Q); % Q is guaranteed to be hermitian since it is a real
%symmetric matrix
% Assume vals is a 1xn matrix
% vecs is an mxn matrix
% both are already sorted
% eig_ind is the eigenvalue/vector to find the perturbed correction for

val2 = vals(eig_ind);
vec2 = vecs(:,eig_ind);

if issparse(Qp)
    vec2 = sparse(vec2);
end

E1 = vec2' * Qp * vec2;
    