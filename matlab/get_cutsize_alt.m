function cutsize = get_cutsize_alt(A,split_ind)
% Get cutsize for a given partition
% A: Adjacency matrix, with columns and rows reordered in the order we used
% for the 1D placement. This reordering simplifies the cutsize calculation
% split_ind: the index of place1d at which to bipartition the design. Since
% A has been reordered in place1d order, we can use split_ind directly to
% split A into quadrants.
% the top left and bottom right represent intra-partition connections in the
% left and right partitions, respectively
% The top right and bottom left quadrants of the adjacency matrix represent
% inter-partition connections, so to calculate cutsize we can just sum up
% all of nonzero values in either quadrant (we only need one quadrant,
% since they both represent the same connections)

% Sum up all entries in the top-right quadrant, which represents
% inter-partition connections
cutsize = sum(sum( A(1:split_ind,split_ind+1:end) ));
cutsize = full(cutsize);