function [Q vals vecs D A time_parse time_eig node_map] = construct_laplacian_from_hgr_blacklist(filename,num_eigs,blacklist)

tic
% read in the netlist
[Q D A] = parse_hgr_sparse_alt3_blacklist(filename,blacklist);
time_parse = toc;
%disp('...HGR Parsing done!')

%% strip out zero rows for faster eigs
indices = 1:length(Q(1,:));
zero_cols = (sum(abs(Q),1) == 0);
keep_cols = ~zero_cols;

% The node map keeps track of which actual node corresponds to the reduced
% node representation. I.E node 5 in the reduced node rep would actually be
% node 10 if nodes 5-9 in the original were in the blackist
node_map = indices(keep_cols);

Q(:,zero_cols) = []; % strip columns
Q(zero_cols,:) = []; % strip rows

qsize = size(Q);
if (num_eigs > qsize(1))
    num_eigs = qsize(1);
end

%% Get the sorted eigenvalues and eigenvectors
tic 
[vals vecs] = get_sorted_eigs(Q,num_eigs);
time_eig = toc;
%disp('...Eigenvalues found!')