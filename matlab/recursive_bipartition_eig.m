function [terminals blocks] = recursive_bipartition_eig(filename,area_constraint)

num_eigs = 10;
node_areas = 1;
% area_constraint = 0.50;

% [FIX] if we just do blacklist we will end up with laplacians
% that have no entries in some rows/columns
% We could just rip out the zero rows/columns, but that will change the
% node indices, so we need to be able to transform from the condensed to
% the full representations

% for each partitioning level
next_partition_level = 1;
num_blocks_next_level = 1; % start by partitioning entire design into two blocks

while( num_blocks_next_level > 0 )
    % Update partition level IDs
    partition_level = next_partition_level;
    dstr = sprintf('Partition level: %d',partition_level');
    sprintf(dstr);
    next_partition_level = partition_level + 1;
    num_blocks = num_blocks_next_level;
    
    num_blocks_next_level = 0;
    % for each block on this partitioning level
    for block_ind = 1:num_blocks
        
        if (partition_level > 1)
            nodes_in_block = blocks{partition_level}{block_ind};
        else
            nodes_in_block = [1 2 3]; % temp
        end

        % We can only partition blocks that have two or more nodes in them
        if ((length(nodes_in_block) > 4) || (partition_level == 1))
            
            % Allow all nodes at the first level
            % At subsequent partitioning levels, figure out which nodes to
            % exclude
            if (partition_level == 1)
                blacklist = [];
            else
                blacklist = construct_blacklist(nodes_in_block,all_nodes);
            end

            [metrics times matrices eigs partitions] = eig_partitioner_blacklist(filename,num_eigs,node_areas,area_constraint,blacklist);

            if(partition_level == 1)
                adjacency_full = matrices.adjacency;
                all_nodes = 1:length(adjacency_full);
                blocks{partition_level}{1} = all_nodes;
                terminals{partition_level}{1} = 0;
            end
            
            % Store off partitioned nodes
            % Assuming here that we're bipartitioning
            blocks{next_partition_level}{2*block_ind-1} = partitions{1};
            blocks{next_partition_level}{2*block_ind} = partitions{2};

            % Get terminal count for each subblock from this partition
            %terminals{next_partition_level} = -1*ones(1,2*num_blocks); % Initialize terminal list
            terminals{next_partition_level}(2*block_ind-1) = get_cutsize_blacklist(partitions{1}, [partitions{2} blacklist], adjacency_full);
            terminals{next_partition_level}(2*block_ind)   = get_cutsize_blacklist(partitions{2}, [partitions{1} blacklist], adjacency_full);

            num_blocks_next_level = num_blocks_next_level + 2;
        end
    end
end