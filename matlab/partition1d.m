function [ratio_cut_min rcm_ind cutsize_min cm_ind ratio_cut_vec cutsize_vec time_partition partition_ratio] = partition1d(place1d,A,area_constraint,node_areas)
% partition a 1d placement with some area constraint
% place1d: 1xn vector representing the 1D placement (obtained after sorting
% the second smallest eigenvector of the laplacian)
% A: The adjacency matrix
% area_constraint: the smallest or largest area allowed for partitioning 
% the design (as a ratio to total size)
% node_areas: either a scalar (in which case we assume all nodes have equal
% area), or a vector, with the ith element being the area of the ith node
% of the system.
% NOTE: unequal node areas not implemented yet. You can pass something in,
% but this won't do anything with it!

tic

%% reorder adjacency matrix to simplify cutsize calculation
% Doing this splits the adjacency matrix into four quadrants
% the top left and bottom right represent intra-partition connections in the
% left and right partitions, respectively
% The top right and bottom left quadrants of the adjacency matrix represent
% inter-partition connections, so to calculate cutsize we can just sum up
% all of nonzero values in either quadrant (we only need one quadrant,
% since they both represent the same connections)
A_reord = A(place1d,place1d);

%% If we pass in a scalar, assume all nodes have equal area
% [FIX] Actually do something with node areas if we have inequal node areas
if(length(node_areas) == 1) % all nodes have equal area
    node_areas = ones(size(place1d));
end
    
%% Figure out what partitions to check
% Use the area constraint to figure out the most unbalanced partitions we
% will consider, and check those and all partitions between them
split_ind = ceil(area_constraint*length(place1d));

split_ind_min = split_ind;
split_ind_max = length(place1d) - split_ind_min;

if(split_ind_max < split_ind_min); % swap them around!
    temp = split_ind_max;
    split_ind_max = split_ind_min;
    split_ind_min = temp;
end

% Normalize split index to system size, so we can figure out the area
% percentage on each side
partition_ratio = (split_ind_min:split_ind_max)/length(place1d);

%% Check all the partitions and figure out which one is the best!
split_indices = split_ind_min:split_ind_max;
cutsize_vec = -1*ones(1,length(split_indices));
ratio_cut_vec = -1*ones(1,length(split_indices));

for split_ind = split_indices
    
    % [FIX] Faster way to do this would be to loop within ratio_cut/cutsize
    % calculation. Pass the adjacency matrix in ONCE and then just schmoo
    % the quadrant size based on the split_indices vector
    ii = split_ind - split_ind_min + 1;
    [ratio_cut cutsize] = get_ratio_cut_alt(A_reord,place1d,split_ind);
    
    ratio_cut_vec(ii) = ratio_cut;
    cutsize_vec(ii) = cutsize;
end

%% Get min cutsize and ratio cut, and their corresponding partition points
[cutsize_min cm_ind] = min(cutsize);
cm_ind = cm_ind + split_ind_min-1;

[ratio_cut_min rcm_ind] = min(ratio_cut);
rcm_ind = rcm_ind + split_ind_min-1;

time_partition = toc;

    