function [Q D A] = parse_hgr_sparse_alt3_blacklist(filename,blacklist)
% parses a netlist for a subset of nodes
% excludes all nodes listed in blacklist

edge_num = 1;

filetext = fileread(filename); % read entire file into array
lines = strread(filetext,'%s','delimiter','\n');
numlines = length(lines);

%% Process first line and allocate space for everything
line_arr = strread(lines{1},'%d','delimiter',' ');
num_nets = line_arr(1);
num_nodes = line_arr(2);

ds = zeros(1,num_nodes);

% We won't know how much space we really need for the adjacency matrix
% until we finish parsing it.
% We could parse through it once just to figure out the total size, or
% we could just try to overallocate and hope for the best (which is what we
% do)
% Here we're assuming that the average net size should be less than 50
% nodes -- if the average size is actually bigger than that, we'll crash
% when we try exceed the bounds of the vector. That would be bad.
ai = zeros(1,num_nodes*50);
aj = zeros(1,num_nodes*50);
as = zeros(1,num_nodes*50);

%% Process main array
for line_ind = 2:numlines % start on line 2, since line 1 is just header info
    line_arr = strread(lines{line_ind},'%d','delimiter',' ');
    
    % remove any elements from the blacklist
    line_arr = line_arr(~ismember(line_arr,blacklist));

    %Check for valid line structure -- should have more than one element
    %and consist of doubles. Just checking first element for simplicity
    is_valid_line = false;
    more_than_one = (length(line_arr) > 1);
    % [FIX] May need to modify this to treat node deletion case
    % since one way to deal with that is to keep nets that were deleted,
    % but with only one dummy node inside
    if(more_than_one) % ignore lines with only one node, since they're not connected to anything
        is_valid_line = true;
    end
    
    if(is_valid_line)
        net_cells = length(line_arr);
        net_weight = 1/(net_cells-1);
        
        % construct diagonal entries - increment the value of each cell's
        % entry in the diagonal value vector
        ds(line_arr) = ds(line_arr) + net_weight*(net_cells-1);
        
        % Prep for adjacency matrix (off-diag) entries
        elchunk_length = net_cells-1;
        
        cells_without_i = zeros(1,elchunk_length);
        
        % Parse through all cells in the net and add all the connections
        % between them to the adjacency vectors
        for i = 1:net_cells %:net_cells-1
            i_cell = line_arr(i);

            cells_without_i(1:end) = [line_arr(1:i-1)' line_arr(i+1:end)'];
            
            edge_vec_start = edge_num;
            edge_vec_end = edge_num + elchunk_length - 1;
            edge_vec = edge_vec_start:edge_vec_end;
            ai(edge_vec) = i_cell;
            aj(edge_vec) = cells_without_i;
            as(edge_vec) = net_weight;
            
            edge_num = edge_num + length(edge_vec);

        end
    end
end

% truncate the vectors to get rid of those extra 0 elements
ai = ai(ai ~= 0);
aj = aj(aj ~= 0);
as = as(1:length(ai));
dij = 1:num_nodes; % since D is diagonal, we just need indices for the diagonals

% generate the sparse adjacency, degree, and laplacian matrices
A = sparse(ai,aj,as,num_nodes,num_nodes);
D = sparse(dij,dij,ds,num_nodes,num_nodes);
Q = D-A;





    