function [terminals blocks cuts] = recursive_bipartition_eig(filename,area_constraint,max_partition_level)

num_eigs = 10;
node_areas = 1;
% area_constraint = 0.50;

if (max_partition_level == -1) % if we pass in -1, just do all the levels
    max_partition_level = 1e5;
end

% [FIX] if we just do blacklist we will end up with laplacians
% that have no entries in some rows/columns
% We could just rip out the zero rows/columns, but that will change the
% node indices, so we need to be able to transform from the condensed to
% the full representations

% for each partitioning level
next_partition_level = 1;
num_blocks_next_level = 1; % start by partitioning entire design into two blocks

while( (num_blocks_next_level > 0) && (next_partition_level < max_partition_level ) )
    % Update partition level IDs
    partition_level = next_partition_level;
    dstr = sprintf('Partition level: %d',partition_level');
    disp(dstr);
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
        if ( ((length(nodes_in_block) > 8) || (partition_level == 1)) && (partition_level < max_partition_level) )
            
            % Allow all nodes at the first level
            % At subsequent partitioning levels, figure out which nodes to
            % exclude
            if (partition_level == 1)
                blacklist = [];
            else
                blacklist = construct_blacklist(nodes_in_block,all_nodes);
            end

            [metrics times matrices eigs partitions] = eig_partitioner_blacklist(filename,num_eigs,node_areas,area_constraint,blacklist);
            
            if(partitions{1}(1) ~= -1) % Don't do anything for bogus partitions
                if(partition_level == 1)
                    adjacency_full = matrices.adjacency;
                    all_nodes = 1:length(adjacency_full);
                    blocks{partition_level}{1} = all_nodes;
                    terminals{partition_level}(1) = 0;
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
    
    % Calculate number of cuts between each set of subblocks
    num_blocks_this_level = length(blocks{partition_level});
    
    if (partition_level == 1)
        cuts{partition_level}(1) = 0;
    else
        cuts{partition_level} = zeros(num_blocks_this_level,num_blocks_this_level);
        for start_block_ind = 1:num_blocks_this_level
            start_block = blocks{partition_level}{start_block_ind};

            for other_block_ind = (start_block_ind+1):num_blocks_this_level
                other_block = blocks{partition_level}{other_block_ind};
                rest_of_nodes = construct_blacklist([start_block other_block],all_nodes);

                cuts{partition_level}(start_block_ind,other_block_ind) = get_cutsize_between_subblocks(start_block, other_block, rest_of_nodes, adjacency_full);
                cuts{partition_level}(other_block_ind,start_block_ind) = cuts{partition_level}(start_block_ind,other_block_ind);
            end
        end
    end
            
end

        
        