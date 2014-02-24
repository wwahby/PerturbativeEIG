function cutsize = get_cutsize_between_subblocks(start_block, other_block, rest_of_nodes, adjacency)

% create reordering vector
vec = [start_block other_block rest_of_nodes];
bound_start = length(start_block);
bound_stop = bound_start + length(other_block);

% Reorder the adjacency matrix so we can simply sum the contributions in
% each subblock to get inter-block cutsize
A = adjacency(vec,vec);

% Sum the adjacency matrix subblock to get connections between start_block
% and other_block
cutsize = full(sum(sum( A(1:bound_start,bound_start+1:bound_stop) )));