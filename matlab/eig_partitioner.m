function [metrics times matrices eigs] = eig_partitioner(filename,num_eigs,node_areas,area_constraint)
% Standard EIG algorithm for spectral partitioning
% filename - character array of the filename containing the .hgr netlist
%    netlist file must be space delimited, no trailing spaces
%    first row must have only two elements: number of nets, and number of
%    nodes, in that order
%    Each subsequent row represents a single net, and has a space-delimited list of connected nodes (indicated by
%    their indices)
% num_eigs - the number of eigenvectors to find from the graph laplacian
% node_areas - 1xn vector with area of each node, OR a scalar (indicating
%    that all nodes have same area). Currently this doesn't actually do
%    anything, and all nodes are assumed to have the same area
% area_constraint - smallest or largest area (as fraction of total area)
%    allowed for one partition
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

%% Parse HGR and construct Q
[Q vals vecs D A time_parse time_eig] = construct_laplacian_from_hgr(filename,num_eigs);

%% 1D Placement
[p1d val2 vec2 time_place] = place_1d(vals,vecs);

%% Partition
[ratio_cut_min rcm_ind cutsize_min cm_ind ratio_cut_vec cutsize_vec time_partition partition_ratio] = partition1d(p1d,A,area_constraint,node_areas);

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
times.eigs = time_eig;
times.total = times.partition + times.placement + times.parse + times.eigs;

% useful matrices
matrices.laplacian = Q;
matrices.degree = D;
matrices.adjacency = A;

%eigenvalues
eigs.vecs = vecs;
eigs.vals = vals;
eigs.val2 = val2;
eigs.vec2 = vec2;

