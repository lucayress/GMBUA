function [M]=match_endmembers(A1, A2)
% Computes the distance between two abundance maps
% A1 and A2: abundance matrices
[P, N] = size(A1);

% compute the pairwise dustances between the endmembers
Cost = zeros(P,P);
for i=1:P
    for j=1:P
        Cost(i,j) = norm(A1(i,:) - A2(j,:), 'fro')/N;
    end
end

% find the cost of the best assignment
% costUnmatched = 0;
%[M,uR,uC] = matchpairs(Cost,costUnmatched,'max');
costUnmatched = -1e6 * P; % large neg value so no endmember is unassigned
[M,uR,uC] = matchpairs(-Cost,costUnmatched,'max');