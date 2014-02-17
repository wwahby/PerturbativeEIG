% Recursive bipartitioning 

% [FIX] if we just do blacklist we will end up with laplacians
% that have no entries in some rows/columns
% We could just rip out the zero rows/columns, but that will change the
% node indices, so we need to be able to transform from the condensed to
% the full representations

% for each partitioning level
next_partition_level = 1; % [FIX] temporary placeholder until we get to outer loop
num_blocks_next_level = 1; % start by partitioning entire design into two blocks

while( num_blocks_next_level > 0 )
    % Update partition level IDs
    partition_level = next_partition_level;
    next_partition_level = partition_level + 1;
    num_blocks = num_blocks_next_level;
    
    % for each block on this partitioning level
    for block_ind = 1:num_blocks

        num_blocks_next_level = 0;

        nodes_in_block = blocks{partition_level}{block_ind};

        % We can only partition blocks that have two or more nodes in them
        if (length(nodes_in_block) > 1)
            blacklist = construct_blacklist(nodes_in_block,all_nodes);

      
            [metrics times matrices eigs partitions] = eig_partitioner_blacklist(filename,num_eigs,node_areas,area_constraint,blacklist);

            % Store off partitioned nodes
            % Assuming here that we're bipartitioning
            blocks{next_partition_level}{2*block_ind-1} = partitions{1};
            blocks{next_partition_level}{2*block_ind} = partitions{2};

            % Get terminal count for each subblock from this partition
            terminals{next_partition_level} = -1*ones(1,2*num_blocks); % Initialize terminal list
            terminals{next_partition_level}(2*block_ind-1) = get_cutsize(nodes_a, adjacency_full);
            terminals{next_partition_level}(2*block_ind)   = get_cutsize(nodes_b, adjacency_full);

            num_blocks_next_level = num_blocks_next_level + 2;
        end
    end
end