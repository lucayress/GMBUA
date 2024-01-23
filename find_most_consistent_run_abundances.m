function [A_bestrun,cost,bestrun_idx]=find_most_consistent_run_abundances(A_all)
% ======================
% Finds the most consistent run
% A_all : cell array containing the abundances, of size P by N, where P 
%         is the number of endmembers in the scene
% 
% A_bestrun : "most representative" abundance from A_all
% ======================

num_runs = length(A_all);

distmtx = zeros(num_runs,num_runs); % adjacency matrix

for i=1:num_runs
    for j=1:num_runs
        if i < j
            distmtx(i,j) = compute_abundance_distance_assignement(A_all{i}, A_all{j});
        else
            distmtx(i,j) = 0; % no self-loops
        end
    end
end
distmtx = distmtx + distmtx';


% compute the most central node/abundance index
% [bestrun_idx]=bestrun_smallestdistance(distmtx);
[bestrun_idx]=bestrun_graphtree(distmtx);

cost = sum(sum(distmtx(bestrun_idx,:)));
A_bestrun = A_all{bestrun_idx};
end





function [ii]=bestrun_smallestdistance(distmtx)
% get the run whose factors are closest to every other run
[~,bestrun_idx] = min(sum(distmtx,2));
ii = bestrun_idx;
end


function [ii] = bestrun_graphtree(distmtx)
% finds minimum spanning tree
G = graph(distmtx);
[T,pred] = minspantree(G);
ucc = centrality(T,'closeness');
[~,ii] = max(ucc);
end


function [mindist]=compute_abundance_distance_assignement(A1, A2)
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

% The cost of the matches in M is sum([Cost(M(1,1),M(1,2)), Cost(M(2,1),M(2,2)), ..., Cost(M(p,1),M(p,2))]).
mindist = 0;
for i=1:size(M,1)
    mindist = mindist + Cost(M(i,1),M(i,2));
end
end




