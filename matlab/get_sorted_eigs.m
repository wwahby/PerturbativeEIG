function [vals vecs] = get_sorted_eigs(Q,num_eigs)
% return 1xn vector of sorted eigenvalues
% and eigenvectors as column vectors in the same order as the eigenvalues

if issparse(Q)
    opts.issym = 1; % The laplacian is symmetric, so use some optimizations
    opts.isreal = 1; % The laplacian is real, so use some optimizations
    opts.disp = 0; % don't spam the console with info about the solution
    [vecs vals] = eigs(Q,num_eigs,-1,opts); % using -1 because in some cases 'SM' or 0 will cause poorly-conditioned matrices during LU factorization
    % This seems to be happening because eigs uses our input as a shift value during factorization (or something)
    % -1 should still be fine since all eigenvalues should be >= 0, so
    % we'll get the same behavior as using 'SM', but without the issues of
    % using 0 as the shift value.
else
    % If the laplacian ISN'T sparse, just use the normal direct solver
    [vecs vals] = eig(Q);
end

vals = diag(vals); % just get nx1 vector of eigenvalues

[vals order] = sort(vals); % get sorted eigenvalues and their order
vecs = vecs(:,order); % get eigenvectors in sorted eigenvalue order