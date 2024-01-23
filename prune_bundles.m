function [S, groups2, bundles] = prune_bundles(S, groups, max_angle)

P = length(unique(groups));
bundles = cell(1,P);
S2 = [];
groups2 = [];
for i = 1:P
    ems = find(groups == i )';
    bundles{1,i} = S(:,ems);
    % plot_bundles(bundles{1,i}, ones(1,size(bundles{1,i},2)));
    [bundles{1,i}, idx] = prune_library_max(bundles{1,i}, max_angle);
    % plot_bundles(bundles{1,i}, ones(1,size(bundles{1,i},2)));
    S2 = [S2, bundles{1,i}];
    groups2 = [groups2; i*ones(size(bundles{1,i},2),1)];
end
S = S2;

end