function [ratio_cut cutsize] = get_ratio_cut_alt(A,place1d,split_ind)
% Get ratio cut for a given partition
% A: Adjacency matrix, with columns and rows reordered in the order we used
% for the 1D placement. This reordering simplifies the cutsize calculation
% place1d: 1xn vector of the 1D placement.
% split_ind: the index of place1d at which to bipartition the design

cutsize = get_cutsize_alt(A,split_ind);
ratio_cut = cutsize/(split_ind * (length(place1d)-split_ind));