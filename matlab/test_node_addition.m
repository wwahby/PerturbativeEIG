%% Test node addition perturbation

filename = '../netlists/ckt_f216.hgr';
filename_p = '../netlists/ckt_f216_nodes_added_1.hgr';

fprintf('\tPartitioning original graph...\n')
[metrics, times, matrices, eigs] = eig_partitioner(filename,num_eigs,node_areas,area_constraint); % normal

fprintf('\tPartitioning perturbed graph (exact)...\n')
[metrics_pe, times_pe, matrices_pe, eigs_pe] = eig_partitioner(filename_p,num_eigs,node_areas,area_constraint); % exact perturbation, only need 2 eigenvalues here


Qo = matrices.laplacian;
Qpe = matrices_pe.laplacian;

Qo_new = zeros(size(Qpe));
Qo_new(1:length(Qo), 1:length(Qo)) = Qo(1:end, 1:end);
Qp = Qpe - Qo_new;

% % Run the perturbed solver
% eigs_to_correct = 2;
% fprintf('\tPartitioning perturbed graph (approximate)...\n')
% [metrics_p, times_p, matrices_p, eigs_p] = eig_partitioner_perturbed(filename_p,matrices.laplacian,eigs.vals,eigs.vecs,node_areas,area_constraint,eigs_to_correct); % approx perturbation
