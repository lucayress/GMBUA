%% Synthetic hyperspectral image generator with endmember P
% --------------------------------------------------------------------
clc, clearvars, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs

load USGS_1995_Library.mat                               % library

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build the dictionary

%  order bands by increasing wavelength
[dummy index] = sort(datalib(:,1));
S =  datalib(index,4:end);
names = names(4:end,:);
L = size(S,1);

% prune the library
min_angle = 4.44;
[S, index] = prune_library2(S,min_angle); % 240  signature
names = names(index',:);

% order  the columns of S by decreasing angles
[S, index, angles] = sort_library_by_angle(S);
names = names(index',:);
namesStr = char(names);

% plot first
n_ems = 25;
groups = [1:n_ems (n_ems+1)*ones(1,size(S,2)-n_ems)];
plot_bundles(S, groups);

ref_ems = [1 3 4 6 7 13 15 18 20 22];
names = names(ref_ems',:);
namesStr = char(names);
ref_ems = S(:,ref_ems);
plot_bundles(ref_ems, 1:size(ref_ems,2));

disp("")


