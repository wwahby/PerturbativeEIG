function [metrics times matrices eigs] = eig_partitioner_perturbed(perturbed_netlist_file,Qorig,vals,vecs,node_areas,area_constraint,eigs_to_correct)
% Runs a modified version of the EIG algorithm
% Uses corrections from TIPT (See Griffiths QM 2nd ed, Ch 6) to calculate
% approximate eigenvectors and eigenvalues for perturbed system
% perturbed_netlist_file - name of file storing perturbed netlist
% Q - original laplacian matrix (sparse)
% vals - (1xnum_vecs) vector of eigenvalues to use for computation
% vecs - (nxnum_vecs) matrix with column vectors of eigenvectors to use for computation
%    vals and vecs must be same length
% node_areas - 1xn vector with area of each node, OR a scalar (indicating
%    that all nodes have same area). Currently this doesn't actually do
%    anything, and all nodes are assumed to have the same area
% area_constraint - smallest or largest area (as fraction of total area)
%    allowed for one partition
% eigs_to_correct - vector of eigenvalue indices to apply approximate perturbation correction to
% ====================================
% Outputs results in container objects
%    metrics
%       cutsize
%          vec - vector of cutsizes for all partitioning solutions
%          min - minimum cutsize for given area constraint
%          min_ind - index of partition cut used for cutsize.min
%          skew.vec - Best cutsize for given area skew
%       ratio_cut
%          vec - vector of ratio cuts for all partitioning solutions
%          min - minimum ratio cuts for given area constraint
%          min_ind - index of partition cut used for ratio_cut.min
%          skew.vec - Best ratio cut for given area skew
%       partition_ratio - size of left partition, as a fraction of total size
%       skew - Area skew corresponding to cutsize.skew.vec and ratio_cut.skew.vec
%   times
%      parse - time to parse netlist file
%      eigs - time to calculate eigenvectors
%      placement - time to perform 1d placement
%      partitioning - time to choose best partitioning solution within area constraint
%      total - total time from start to finish
%   matrices
%      laplacian - the laplacian matrix corresponding to the input netlist
%      adjacency - the adjacency matrix
%      degree - the degree matrix
%   eigs
%      vals - the first num_eigs eigenvalues of the laplacian
%      vecs - the first num_eigs eigenvectors of the laplacian
%      val2 - the second eigenvalue of the laplacian
%      vec2 - the second eigenvector of the laplacian

%% Parsing
% Parse the new netlist file, generate the perturbation to the laplacian
% (for partitioning) and the new adjacency matrix (which will be used for cutsize evaluation)
tic
[Qpe Dpe Ape] = parse_hgr_sparse_alt3(perturbed_netlist_file);
Qpe = (Qpe + Qpe')/2; % Ensure symmetry despite weird rounding errors
time_parse = toc;

%% Determine whether nodes have been added
if (length(Qpe) > length(Qorig))
    num_eigs = length(vals);
    Qb = Qpe(length(Qorig)+1:end, length(Qorig)+1:end); % construct new node laplacian
    
    num_eigs = min(num_eigs, length(Qb) );
    [vals_b, vecs_b] = get_sorted_eigs(Qb,num_eigs);
    %vals_b = diag(vals_b)';
    
    vecs_b_full = sparse(length(Qpe), length(vals_b));
    for vind = 1:length(vecs_b(1,:))
        vecs_b_full(length(Qorig)+1:end, vind) = vecs_b(:,vind);
    end
    
    vals_ab = [vals', vals_b'];
    
    [vals_ab_sorted, vals_ab_sorted_inds] = sort(vals_ab);
    
    vecs_ab = sparse(length(Qpe), length(Qpe));
    vecs_ab(1:end-length(vals_b),1:length(vals)) = vecs;
    vecs_ab(:,length(Qorig)+1:end) = vecs_b_full(:,1:end);
    vecs_ab_sorted = vecs_ab(:, vals_ab_sorted_inds);
    
    Qo_new = sparse(length(Qpe), length(Qpe));
    Qo_new(1:length(Qorig), 1:length(Qorig)) = Qorig(1:end, 1:end);
    
    Qp = Qpe - Qo_new;
    vals = vals_ab_sorted;
    vecs = vecs_ab_sorted;
else
    Qp = Qpe - Qorig;
end

%% 1D Placement
num_vecs_to_use = length(vecs(1,:)); % Use all the eigenvectors we're supplied
[p1d vals_p vecs_p time_place] = place_1d_perturbed_alt(Qp,vals,vecs,num_vecs_to_use,eigs_to_correct);

%% Partition
[ratio_cut_min rcm_ind cutsize_min cm_ind ratio_cut_vec cutsize_vec time_partition partition_ratio] = partition1d(p1d,Ape,area_constraint,node_areas);

%% Calculate skew
[c_min skew] = eval_with_skew(cutsize_vec,partition_ratio);
[rc_min skew] = eval_with_skew(ratio_cut_vec,partition_ratio);

%% Gather outputs

% Metrics (ratio cut, cutsize, and other useful data)
metrics.ratio_cut.vec = ratio_cut_vec;
metrics.ratio_cut.min = ratio_cut_min;
metrics.ratio_cut.min_ind = rcm_ind;
metrics.ratio_cut.skew.vec = rc_min;
metrics.cutsize.vec = cutsize_vec;
metrics.cutsize.min = cutsize_min;
metrics.cutsize.min_ind = cm_ind;
metrics.cutsize.skew.vec = c_min;
metrics.partition_ratio = partition_ratio;
metrics.skew = skew;

% timing data
times.partition = time_partition;
times.placement = time_place;
times.parse = time_parse;
times.eigs = 0; % no eigenvalues calculated
times.total = times.partition + times.placement + times.parse + times.eigs;

% useful matrices
matrices.laplacian = Qp;
matrices.degree = 0; % Degree matrix not available
matrices.adjacency = Ape;

%eigenvalues
eigs.vecs = vecs_p;
eigs.vals = vals_p;
eigs.val2 = vals_p(2); % perturbed 2nd eigenvalue
eigs.vec2 = vecs_p(:,2); % perturbed 2nd eigenvector

