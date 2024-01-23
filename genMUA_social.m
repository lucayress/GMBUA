function [X_hat_l1_spreg] = genMUA_social(Xim,datalib,groups,A_init,slic_size,slic_reg,beta,lambda,rho,maxiter_ADMM,type,fraction,tol_a,rand_state,vlfeat_path,verbose)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define random states
if rand_state ~= 0
    rand('state',rand_state);
    randn('state',rand_state);
end

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

% abundance initialization
A_init_sppx = FCLSU(squeeze(avg_superpx)',A)';

% unmixing
[X_hat_l1_suppx, optim_struct_fractional] = social_unmixing(squeeze(avg_superpx)',A,groups,A_init_sppx,lambda,rho,maxiter_ADMM,type,fraction,tol_a,verbose);

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

% Yc = A*X_hat_l1_spreg;
% Yc_im = reshape(Yc', nl, nc, L);
% Xc_global = bundle2global(X_hat_l1_spreg, A, groups);
% Xc_global_im = reshape(Xc_global', nl, nc, 5);

%% Generalized MUA
Y_new = [Y' (sqrt(beta)*X_hat_l1_spreg)']';
A_new = [A' (sqrt(beta)*eye(size(A,2)))']';

[X_hat_l1_spreg, optim_struct_fractional] = social_unmixing(Y_new,A_new,groups,A_init,lambda,rho,maxiter_ADMM,type,fraction,tol_a,verbose);
