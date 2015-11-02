function [vals vecs] = get_sorted_eigs(Q,num_eigs)
% return 1xn vector of sorted eigenvalues
% and eigenvectors as column vectors in the same order as the eigenvalues
% May not return num_eigs values, as some may be spurious negative eigs

if issparse(Q)
    opts.issym = true; % The laplacian is symmetric, so use some optimizations
    opts.isreal = true; % The laplacian is real, so use some optimizations
    opts.disp = 0; % don't spam the console with info about the solution
%     opts.tol = 1e-20; % really fine tolerance
    opts.maxit = 3e4; % lots of iterations

    Q = (Q + Q')/2; % Ensure symmetry
    [vecs vals] = eigs(Q,num_eigs,'SA',opts); % using -1 because in some cases 'SM' or 0 will cause poorly-conditioned matrices during LU factorization
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

%% throw away negative eigenvalues 
% MATLAB may give us spurious negative eigenvalues -- use this to throw
% them away
% Best method would be to use 'SM' and shift-invert mode to find actual smallest
% eigenvalues, but it occasionally fails because 0 is an eigenvalue
% posvals = (vals > 0);
% vals = vals(posvals);
% vecs = vecs(:,posvals);

%% Throw away really tiny positive values, as they're probably spurious too
bigvals = (vals > 1e-5);

% we need to keep the first "zero" eigenvalue
if(isempty(bigvals)) % do nothing if no good eigenvalues
    vals = -1;
    vecs = -1;
else
    if(bigvals(1) == 0)
        ind = find(bigvals == 0,1,'last');
        bigvals(ind) = 1;
    end
    
    vals = vals(bigvals);
    vecs = vecs(:,bigvals);
    
end


