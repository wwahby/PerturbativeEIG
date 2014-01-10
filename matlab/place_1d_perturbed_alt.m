function [p1dp Ep psip time_place_p] = place_1d_perturbed_alt(Qp,vals,vecs,num_vecs_to_use)
%[vecs vals] = eig(Q); % Q is guaranteed to be hermitian since it is a real
%symmetric matrix
% Assume vals is a 1xn matrix
% vecs is an mxn matrix
% both are already sorted

tic
E1 = calc_perturbed_eigval_first_order(Qp,vals,vecs);
[psi1_full psi1] = calc_perturbed_eigvec_first_order(Qp,vals,vecs,num_vecs_to_use);
E2 = calc_perturbed_eigval_second_order(Qp,vals,vecs,num_vecs_to_use);
psi2 = calc_perturbed_eigvec_second_order(Qp,vals,vecs,num_vecs_to_use,E1,psi1);

psip = vecs(:,2) + psi1 + psi2;
% renormalize psip
%psip = psip./(psip'*psip);

Ep = vals(2) + E1 + E2;

[psip_sorted p1dp] = sort(psip);
time_place_p = toc;