function cutsize = get_cutsize_blacklist(block,other_nodes,adjacency)

vec = [block other_nodes];
bound = length(block);
A = adjacency(vec,vec);

cutsize = full(sum(sum( A(1:bound,bound+1:end) )));