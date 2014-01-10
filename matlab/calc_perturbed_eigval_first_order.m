function [E1 ] = calc_perturbed_eigval_first_order(Qp,vals,vecs)
%[vecs vals] = eig(Q); % Q is guaranteed to be hermitian since it is a real
%symmetric matrix
% Assume vals is a 1xn matrix
% vecs is an mxn matrix
% both are already sorted

val2 = vals(2);
vec2 = vecs(:,2);

if issparse(Qp)
    vec2 = sparse(vec2);
end

E1 = vec2' * Qp * vec2;
    