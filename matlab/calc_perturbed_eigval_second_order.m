function [E2] = calc_perturbed_eigval_second_order(Qp,vals,vecs,num_vecs_to_use,eig_ind)
%[vecs vals] = eig(Q); % Q is guaranteed to be hermitian since it is a real
%symmetric matrix
% Assume vals is a 1xn matrix
% vecs is an mxn matrix
% both are already sorted

val2 = vals(eig_ind);
vec2 = vecs(:,eig_ind);

% if issparse(Qp)
%     vec2 = sparse(vec2);
% end

E2 = 0;

if(num_vecs_to_use < 0)
    num_vecs_to_use = length(vals);
end

vecs_to_include = [1:(eig_ind-1) (eig_ind+1):num_vecs_to_use];

top1 = Qp*vec2;
for i = vecs_to_include
    vv = vecs(:,i);
    top = vv' * top1;
    
    bot = val2 - vals(i);
    
    E2 = E2 + abs(top)^2/bot;
end
    
    