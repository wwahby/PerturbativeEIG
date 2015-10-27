function [Q, vals, vecs, D, A, time_parse, time_eig] = construct_laplacian_from_hgr(filename,num_eigs)

tic
% read in the netlist
[Q, D, A] = parse_hgr_sparse_alt3(filename);
time_parse = toc;
%disp('...HGR Parsing done!')

sufficient_vals = 0;
num_runs = 0;
max_runs = 10;

while ((sufficient_vals == 0) && (num_runs < max_runs))
    num_runs = num_runs + 1;
    % Get the sorted eigenvalues and eigenvectors
    tic
    [vals, vecs] = get_sorted_eigs(Q,num_eigs);
    time_eig = toc;

    if (length(vals) > 2)
        sufficient_vals = 1;
    else
        fprintf('\t\tNot enough true eigenvalues found! Rerunning (%d/%d)...\n',num_runs, max_runs);
    end
end