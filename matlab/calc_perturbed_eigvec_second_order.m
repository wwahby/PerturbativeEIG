function psi2 = calc_perturbed_eigvec_second_order(Qp,vals,vecs,num_vecs_to_use,E1,psi1)
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

val_prt = 0;

if(num_vecs_to_use < 0)
    num_vecs_to_use = length(vals);
end

vecs_to_include = [1 3:num_vecs_to_use];

for i = vecs_to_include
    psii = vecs(:,i);
    c(i) = (E1* (psii'*psi1)  - (psii'*(Qp*psi1) ));
    c(i) = c(i)/(vals(i) - E1);
end

psi2 = zeros(size(vec2));
for i = vecs_to_include
    psi2 = psi2 + c(i)*vecs(:,i);
end

%psi2_full = vec2 + psi2;
    
    