function [vec_min skew] = eval_with_skew(vec,partition_ratio)
% Convert a full range of partitioning solutions into a list of the best
% partitioning options for each value of area skew.
%    i.e. a 60/40 split and a 40/60 split of the original design are both
%    partitions with 10% area skew, so choose the one that gives the best
%    result
% Takes in a vector of cutsizes or ratio-cuts, ordered from smallest
%    left-partition size to largest left-partition size (center point should
%    be a balanced partition)
% Finds the best partitioning solution for each area constraint

vec_min = min(vec,fliplr(vec));
mid = round(length(vec_min)/2);
vec_min = vec_min(mid:end);
skew = partition_ratio(mid:end)-0.5;