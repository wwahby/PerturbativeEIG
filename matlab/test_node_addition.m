%% Test node addition perturbation

filename = '../netlists/ckt_f216.hgr';
filename_p = '../netlists/ckt_f216_nodes_added_1.hgr';
filename_p = '../netlists/ckt_f216_nodes_added_3.hgr';

num_eigs = 7;
node_areas = 1;
area_constraint = 0.45;

fprintf('\tPartitioning original graph...\n')
[metrics, times, matrices, eigs] = eig_partitioner(filename,num_eigs,node_areas,area_constraint); % normal

fprintf('\tPartitioning perturbed graph (exact)...\n')
[metrics_pe, times_pe, matrices_pe, eigs_pe] = eig_partitioner(filename_p,num_eigs,node_areas,area_constraint); % exact perturbation, only need 2 eigenvalues here


Qo = matrices.laplacian;
Qpe = matrices_pe.laplacian;

[vecs, vals] = eig(full(Qo));
vals = diag(vals)';
[vecs_pe, vals_pe] = eig(full(Qpe));
vals_pe = diag(vals_pe)';

Qb = full(Qpe(length(Qo)+1:end, length(Qo)+1:end));

[vecs_b, vals_b] = eig(Qb);
vals_b = diag(vals_b)';
vecs_b_full = zeros(length(Qpe), length(vals_b));

for vind = 1:length(vecs_b(1,:))
    vecs_b_full(length(Qo)+1:end, vind) = vecs_b(:,vind);
end

vals_ab = [vals, vals_b];

[vals_ab_sorted, vals_ab_sorted_inds] = sort(vals_ab);

vecs_ab = zeros(length(Qpe), length(Qpe));
vecs_ab(1:end-length(vals_b),1:length(vals)) = vecs;
vecs_ab(:,length(Qo)+1:end) = vecs_b_full(:,1:end);
vecs_ab_sorted = vecs_ab(:, vals_ab_sorted_inds);




Qo_new = zeros(size(Qpe));
Qo_new(1:length(Qo), 1:length(Qo)) = Qo(1:end, 1:end);
Qp = Qpe - Qo_new;

eigs_to_correct = 2;
num_vecs_to_use = 7;
[vecs_p, vals_p] = approximate_perturbed_eigs(Qp,vecs_ab_sorted,vals_ab_sorted,eigs_to_correct,num_vecs_to_use);


err_vals_rel = abs(vals_p - vals_pe)./vals_pe;
[vals_sorted_p, order_p] = sort(vecs_p(:,2));
[vals_sorted_pe, order_pe] = sort(vecs_pe(:,2));


% % Run the perturbed solver
% eigs_to_correct = 2;
% fprintf('\tPartitioning perturbed graph (approximate)...\n')
% [metrics_p, times_p, matrices_p, eigs_p] = eig_partitioner_perturbed(filename_p,matrices.laplacian,eigs.vals,eigs.vecs,node_areas,area_constraint,eigs_to_correct); % approx perturbation
