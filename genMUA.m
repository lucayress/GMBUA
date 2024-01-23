function [X_hat_l1_spreg] = genMUA(Xim, datalib, slic_size, slic_reg, lambda1_sp, lambda2_sp, beta, maxiter_ADMM, tol_a, rand_state, vlfeat_path)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define random states
rand('state',rand_state);
randn('state',rand_state);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load fractional abundances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Size of the images
nl = size(Xim,1);
nc = size(Xim,2);

%% Dictionary

A =  datalib;
L = size(A,1);
n = size(A,2);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiscale regularization using superpixels SLIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Decompose the image into homogeneous regions and apply clustering
% individually

% Set VL_feat toolbox
run(strcat(vlfeat_path,'vlfeat-0.9.21/toolbox/vl_setup'));

X = reshape(Xim,nl*nc,L)';
Y = X;
Y2 = reshape(X', nl, nc, L);
Y2a = Y2;

% tic

% reorder and rescale data into 2-D array
[numRows,numCols,numSpectra] = size(Y2);
scfact = mean(reshape(sqrt(sum(Y2.^2,3)), numRows*numCols, 1));
Y2 = Y2./scfact;

% compute superpixels
spSegs = vl_slic(single(Y2), slic_size, slic_reg);
numSuperpixels = double(max(spSegs(:)))+1;
% ------
% Unmix the superpixels

avg_superpx = zeros(1, numSuperpixels+1, L);

for i=0:numSuperpixels
    [rowi, coli] = find(spSegs==i);
    
    for j=1:length(rowi)
        % Averages all pixels inside each superpixel
        if j == 1
            avg_superpx(1,i+1,:) = (1/length(rowi)) * Y2a(rowi(j),coli(j),:);
        else
            avg_superpx(1,i+1,:) = avg_superpx(1,i+1,:) + (1/length(rowi)) * Y2a(rowi(j),coli(j),:);
        end
    end
end

% Unmix each superpixel individually
[X_hat_l1_suppx] = sunsal(A,squeeze(avg_superpx)','lambda',lambda1_sp,'ADDONE','no','POSITIVITY','yes', ...
    'TOL',tol_a, 'AL_iters',maxiter_ADMM,'verbose','off');

% Re-attribute the abundances for the entire matrix
temp = zeros(size(Y2,1), size(Y2,2), n);
for i=0:numSuperpixels
    [rowi, coli] = find(spSegs==i);
    
    % Attributes unmixing result to all pixels in a voxel
    for j=1:length(rowi)
        temp(rowi(j),coli(j),:) = X_hat_l1_suppx(:,i+1);
    end
end

X_hat_l1_spreg = reshape(temp, [size(Y2,1)*size(Y2,2) n])';

%% Generalized MUA
Y_new = [Y' (sqrt(beta)*X_hat_l1_spreg)']';
A_new = [A' (sqrt(beta)*eye(size(A,2)))']';
 
[X_hat_l1_spreg] = sunsal(A_new, Y_new,'lambda',lambda2_sp,'ADDONE','no','POSITIVITY','yes', ...
                       'TOL',tol_a, 'AL_iters',maxiter_ADMM,'verbose','no');
