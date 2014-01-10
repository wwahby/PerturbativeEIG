function z = calc_l_squared(A,vec)
% z = calc_l_squared(A,vec)
% A is the adjacency matrix for a graph representing a netlist
% vec is the 2nd-smallest eigenvector (unsorted) of the corresponding
% Laplacian matrix
% This function calculates the weighted l^2 of the elements in vec
% z = 1/2 sum_i,j { a_ij * l_ij^2 }
% the 1/2 in front accounts for the fact that A is symmetric so a_ij = a_ji
% and l_ij = l_ji by definition

% We can automatically incorporate the 1/2 by only iterating over half the
% matrix
% We can also avoid the diagonal elements since A has only zeros on the
% diagonal
% 
% [m n] = size(A);
% 
% z = 0;
% for i = 1:m
%     for j = i+1:n
%         lij = vec(j) - vec(i);
%         z = z + A(i,j)*lij^2;
%     end
% end


z = 0;
[i j a] = find(A); % get indices and values of elements in A

for k = 1:length(a)
    lij = vec(j(k)) - vec(i(k));
    z = z + a(k)*lij^2;
end
z = 1/2*z;