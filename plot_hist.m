function plot_hist(groups)
% Plot histogram of no. of signatures per bundle

P = length(unique(groups));
edges = 0.5:1:P+0.5;
figure; set(gcf,'color', 'white')
histogram(groups,edges), xticks(1:P)

end