function plot_bundles(bundle, groups, M_sort)

S = bundle;
L = size(S,1);
P = length(unique(groups));

if nargin == 2
    M_sort = (1:P)';
end

if P == 1
    cols = 1;
elseif P < 5
    cols = P;
else
    cols = 5;
end
rows = ceil(P/5);
vertical_spacing = 0.001;
horizontal_spacing = 0.02;
figure; set(gcf, 'Color', 'w');
for i=1:P
    ems = find( groups == M_sort(i,1) );
    subaxis(rows,cols,i,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    for j=1:length(ems)
        plot(1:L,S(:,ems(j)),'linewidth',1)
        hold on
    end
    title("Bundle " + i)
    xlim([1 L]); ylim([0 1.05]);
    grid on; axis square;
end

end