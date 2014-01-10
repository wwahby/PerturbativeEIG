function [p1d val2 vec2 time_place] = place_1d(vals,vecs)
% Takes sorted list of eigenvalues and eigenvectors (each column of vec
% represents an eigenvector) and returns a 1D placement (p1d) of nodes involved
% based on the second smallest eigenvector
% Assume vals and vecs already sorted


tic
val2 = vals(2);
vec2 = vecs(:,2);

% we just care about the sorting order -- we don't really need the
% eigenvector after it's been sorted, so just throw it away
[vec_sorted p1d] = sort(vec2);

time_place = toc;