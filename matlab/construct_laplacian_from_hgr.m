function [Q vals vecs D A time_parse time_eig] = construct_laplacian_from_hgr(filename,num_eigs)

tic
% read in the netlist
[Q D A] = parse_hgr_sparse_alt3(filename);
time_parse = toc;
%disp('...HGR Parsing done!')

tic
% Get the sorted eigenvalues and eigenvectors
[vals vecs] = get_sorted_eigs(Q,num_eigs);
time_eig = toc;
%disp('...Eigenvalues found!')