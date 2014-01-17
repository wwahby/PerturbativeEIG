function [p1dp vals_p vecs_p time_place_p] = place_1d_perturbed_alt(Qp,vals,vecs,num_vecs_to_use,eigs_to_correct)
%[vecs vals] = eig(Q); % Q is guaranteed to be hermitian since it is a real
%symmetric matrix
% Assume vals is a 1xn matrix
% vecs is an mxn matrix
% both are already sorted

tic
%eigs_to_correct = 2; % correct 2nd smallest eigvec/val
% E1 = calc_perturbed_eigval_first_order(Qp,vals,vecs,eig_ind);
% [psi1_full psi1] = calc_perturbed_eigvec_first_order(Qp,vals,vecs,num_vecs_to_use,eig_ind);
% E2 = calc_perturbed_eigval_second_order(Qp,vals,vecs,num_vecs_to_use,eig_ind);
% psi2 = calc_perturbed_eigvec_second_order(Qp,vals,vecs,num_vecs_to_use,E1,psi1,eig_ind);
[vecs_p vals_p] = approximate_perturbed_eigs(Qp,vecs,vals,eigs_to_correct,num_vecs_to_use);

psip = vecs_p(:,2);
% renormalize psip
%psip = psip./(psip'*psip);

Ep = vals_p(2);

[psip_sorted p1dp] = sort(psip);
time_place_p = toc;