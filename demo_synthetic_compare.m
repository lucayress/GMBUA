% This script reproduces the results of the algorithms presented in
%
%  A Generalized Multiscale Bundle-Based Hyperspectral Sparse Unmixing Algorithm
%  L. C. Ayres, R.A. Borsoi, J.C.M. Bermudez, S.J.M. de Almeida.
%  IEEE Geoscience and Remote Sensing Letters, 2023.
%
%
% Author: Luciano Ayres
% Latest Revision: 23-Jan-2024
% Revision: 1.0

%% DEMO - GMBUA - Synthetic data

clc, clearvars, close all

%% Load data

path = 'test_data/';
image_name = 'synthetic';
vlfeat_path = '';

load(strcat('test_data/',image_name,'/',image_name,'.mat')) % observed data
bundle_true = bundle;
groups_true = groups;

% [A_true_final,S_true] = bundle2global(A_new,bundle,groups);
A_true_final = A_new;
H_true = X;

% Get observed data dimensions
data = Xim;
bundles = 5;
[m,n,L] = size(data);
N = m*n;
X = reshape(data,N,L)';

% Add noise
noise = 15;
if noise ~= 0 % noise in dB 
    sigma = sqrt(sum(sum((X).^2))/N/L/10^(noise/10));
    noise_ = sigma*randn(L,N);
    X = X + noise_;
    Xim = reshape(X', m, n, L);
end

% Show image
figure, set(gcf,'color', 'white')
imshow(Xim(:,:,[29 15 12]))

%% Parameters

% Algorithm runs
K  = 30;    % GMBUA
R  = 30;    % Select the most representative solutions from GMBUA runs

K_ = R;     % others algorithms

% endmember extraction
P = bundles;            % number of endmembers
bundle_nbr = 10;        % number of VCA runs
percent = 10;           % percentage of pixels considered in each run

% ADMM optimization
rho = 10;
tol_a = 10^(-6);
maxiter_ADMM = 1000;
verbose = 0;

%% Variable initialization

% endmember extraction
groups = cell(1,K_);
bundle = cell(1,K_);

% all K abundance estimations
A_FCLSU_all = cell(1,K_);
A_FCLSU_global_all = cell(1,K_);
A_collaborative_all = cell(1,K_);
A_collaborative_global_all = cell(1,K_);
A_group_all = cell(1,K_);
A_group_global_all = cell(1,K_);
A_elitist_all = cell(1,K_);
A_elitist_global_all = cell(1,K_);
A_fractional_all = cell(1,K_);
A_fractional_global_all = cell(1,K_);

M_FCLSU_k = cell(1,K_);
M_collaborative_k = cell(1,K_);
M_group_k = cell(1,K_);
M_elitist_k = cell(1,K_);
M_fractional_k = cell(1,K_);

time_FCLSU = zeros(1,K_);
time_collaborative = zeros(1,K_);
time_group = zeros(1,K_);
time_elitist = zeros(1,K_);
time_fractional = zeros(1,K_);

fprintf('\nFCLSU, collaborative, group, elitist, fractional')
fprintf('\nK = ')
for k = 1:K_
    rng(k);
    fprintf([num2str(k) ', '])
    [groups{1,k}, bundle{1,k}] = batchvca2(X, P, bundle_nbr, percent);   % extract endmember bundles
    % bundle{1,k} = bundle_true;
    % groups{1,k} = groups_true;

    % FCLSU
    tic
    A_FCLSU_all{1,k} = FCLSU(X, bundle{1,k})';
    time_FCLSU(1,k) = toc;
    [A_FCLSU_global_all{1,k}, S_FCLSU_global] = bundle2global(A_FCLSU_all{1,k},bundle{1,k},groups{1,k});
    A_init = A_FCLSU_all{1,k};
    M_FCLSU_k{1,k} = match_endmembers(A_FCLSU_global_all{1,k}, A_true_final);
    A_FCLSU_global_all{1,k} = A_FCLSU_global_all{1,k}(M_FCLSU_k{1,k}(:,1),:);

    % Collaborative sparsity
    type = 'asc';
    lambda = 0.01;
    tic
    [A_collaborative_all{1,k}] = ADMM_collaborative_unmixing(X,A_init,bundle{1,k},lambda,rho,maxiter_ADMM,type,tol_a,verbose);
    time_collaborative(1,k) = toc;
    [A_collaborative_global_all{1,k}, S_collaborative_global] = bundle2global(A_collaborative_all{1,k},bundle{1,k},groups{1,k});
    M_collaborative_k{1,k} = match_endmembers(A_collaborative_global_all{1,k}, A_true_final);
    A_collaborative_global_all{1,k} = A_collaborative_global_all{1,k}(M_collaborative_k{1,k}(:,1),:);

    % Group penalty
    lambda = 0.1;
    fraction = 0.01;
    type = 'group';
    tic
    [A_group_all{1,k}, ~] = social_unmixing(X,bundle{1,k},groups{1,k},A_init,lambda,rho,maxiter_ADMM,type,fraction,tol_a,verbose);
    time_group(1,k) = toc;
    [A_group_global_all{1,k}, S_group_global] = bundle2global(A_group_all{1,k},bundle{1,k},groups{1,k});
    M_group_k{1,k} = match_endmembers(A_group_global_all{1,k}, A_true_final);
    A_group_global_all{1,k} = A_group_global_all{1,k}(M_group_k{1,k}(:,1),:);

    % Elitist penalty
    lambda = 0.1;
    fraction = 0.01;
    type = 'elitist';
    tic
    [A_elitist_all{1,k}, ~] = social_unmixing(X,bundle{1,k},groups{1,k},A_init,lambda,rho,maxiter_ADMM,type,fraction,tol_a,verbose);
    time_elitist(1,k) = toc;
    [A_elitist_global_all{1,k}, S_elitist_global] = bundle2global(A_elitist_all{1,k},bundle{1,k},groups{1,k});
    M_elitist_k{1,k} = match_endmembers(A_elitist_global_all{1,k}, A_true_final);
    A_elitist_global_all{1,k} = A_elitist_global_all{1,k}(M_elitist_k{1,k}(:,1),:);

    % Fractional penalty
    lambda = 0.01;
    fraction = 1;
    type = 'fractional';
    tic
    [A_fractional_all{1,k}, ~] = social_unmixing(X,bundle{1,k},groups{1,k},A_init,lambda,rho,maxiter_ADMM,type,fraction,tol_a,verbose);
    time_fractional(1,k) = toc;
    [A_fractional_global_all{1,k}, S_fractional_global] = bundle2global(A_fractional_all{1,k},bundle{1,k},groups{1,k});
    M_fractional_k{1,k} = match_endmembers(A_fractional_global_all{1,k}, A_true_final);
    A_fractional_global_all{1,k} = A_fractional_global_all{1,k}(M_fractional_k{1,k}(:,1),:);
end
% Find most representative abundance
[A_FCLSU_final, ~, FCLSU_bestrun_idx] = find_most_consistent_run_abundances(A_FCLSU_global_all(1,:));    
A_FCLSU = A_FCLSU_all{1, FCLSU_bestrun_idx};
[A_collaborative_final, ~, collaborative_bestrun_idx] = find_most_consistent_run_abundances(A_collaborative_global_all(1,:));    
A_collaborative = A_collaborative_all{1, collaborative_bestrun_idx};
[A_group_final, ~, group_bestrun_idx] = find_most_consistent_run_abundances(A_group_global_all(1,:));    
A_group = A_group_all{1, group_bestrun_idx};
[A_elitist_final, ~, elitist_bestrun_idx] = find_most_consistent_run_abundances(A_elitist_global_all(1,:));    
A_elitist = A_elitist_all{1, elitist_bestrun_idx};
[A_fractional_final, ~, fractional_bestrun_idx] = find_most_consistent_run_abundances(A_fractional_global_all(1,:));    
A_fractional = A_fractional_all{1, fractional_bestrun_idx};

% Compute metrics

% Reconstruction observed data
H_FCLSU = bundle{1,FCLSU_bestrun_idx}*A_FCLSU;
H_collaborative = bundle{1,collaborative_bestrun_idx}*A_collaborative;
H_group = bundle{1,group_bestrun_idx}*A_group; 
H_elitist = bundle{1,elitist_bestrun_idx}*A_elitist; 
H_fractional = bundle{1,fractional_bestrun_idx}*A_fractional; 

% RMSE
RMSE_FCLSU = sqrt(1/L*sum((H_FCLSU-X).^2,1));
RMSE_collaborative =  sqrt(1/L*sum((H_collaborative-X).^2,1));
RMSE_group = sqrt(1/L*sum((H_group-X).^2,1));
RMSE_elitist = sqrt(1/L*sum((H_elitist-X).^2,1));
RMSE_fractional = sqrt(1/L*sum((H_fractional-X).^2,1));

% SRE
SRE_FCLSU_k = zeros(1,K_);
SRE_collaborative_k = zeros(1,K_);
SRE_group_k = zeros(1,K_);
SRE_elitist_k = zeros(1,K_);
SRE_fractional_k = zeros(1,K_);
for k=1:K_
    SRE_FCLSU_k(1,k) = 20*log10(norm(A_true_final,'fro')/norm(A_FCLSU_global_all{1,k} - A_true_final,'fro'));
    SRE_collaborative_k(1,k) =  20*log10(norm(A_true_final,'fro')/norm(A_collaborative_global_all{1,k} - A_true_final,'fro'));
    SRE_group_k(1,k) = 20*log10(norm(A_true_final,'fro')/norm(A_group_global_all{1,k} - A_true_final,'fro'));
    SRE_elitist_k(1,k) = 20*log10(norm(A_true_final,'fro')/norm(A_elitist_global_all{1,k} - A_true_final,'fro'));
    SRE_fractional_k(1,k) = 20*log10(norm(A_true_final,'fro')/norm(A_fractional_global_all{1,k} - A_true_final,'fro'));
end
SRE_FCLSU = SRE_FCLSU_k(1,FCLSU_bestrun_idx);
SRE_collaborative =  SRE_collaborative_k(1,collaborative_bestrun_idx);
SRE_group = SRE_group_k(1,group_bestrun_idx);
SRE_elitist = SRE_elitist_k(1,elitist_bestrun_idx);
SRE_fractional = SRE_fractional_k(1,fractional_bestrun_idx);

aux = sort(SRE_FCLSU_k);
SRE_FCLSU_median = aux(ceil(length(SRE_FCLSU_k)/2));
FCLSU_median_idx = find(SRE_FCLSU_k == SRE_FCLSU_median, 1, 'first');

aux = sort(SRE_collaborative_k);
SRE_collaborative_median = aux(ceil(length(SRE_collaborative_k)/2));
collaborative_median_idx = find(SRE_collaborative_k == SRE_collaborative_median, 1, 'first');

aux = sort(SRE_group_k);
SRE_group_median = aux(ceil(length(SRE_group_k)/2));
group_median_idx = find(SRE_group_k == SRE_group_median, 1, 'first');

aux = sort(SRE_elitist_k);
SRE_elitist_median = aux(ceil(length(SRE_elitist_k)/2));
elitist_median_idx = find(SRE_elitist_k == SRE_elitist_median, 1, 'first');

aux = sort(SRE_fractional_k);
SRE_fractional_median = aux(ceil(length(SRE_fractional_k)/2));
fractional_median_idx = find(SRE_fractional_k == SRE_fractional_median, 1, 'first');

% SAM
SAM_true = zeros(N,1);
SAM_FCLSU = zeros(N,1);
SAM_collaborative = zeros(N,1);
SAM_group = zeros(N,1);
SAM_elitist = zeros(N,1);
SAM_fractional = zeros(N,1);

for k = 1:N
    SAM_FCLSU(k) = 180/pi*real(acos((X(:,k)'*H_FCLSU(:,k))...
        /(norm(X(:,k))*norm(H_FCLSU(:,k)))));
end

for k = 1:N
    SAM_collaborative(k) = 180/pi*real(acos((X(:,k)'*H_collaborative(:,k))...
        /(norm(X(:,k))*norm(H_fractional(:,k)))));
end

for k = 1:N
    SAM_group(k) = 180/pi*real(acos((X(:,k)'*H_group(:,k))...
        /(norm(X(:,k))*norm(H_group(:,k)))));
end

for k = 1:N
    SAM_elitist(k) = 180/pi*real(acos((X(:,k)'*H_elitist(:,k))...
        /(norm(X(:,k))*norm(H_elitist(:,k)))));
end

for k = 1:N
    SAM_fractional(k) = 180/pi*real(acos((X(:,k)'*H_fractional(:,k))...
        /(norm(X(:,k))*norm(H_fractional(:,k)))));
end

%% Generalized MUA

groups = cell(R,K);
bundle = cell(R,K);
A_MUA_all = cell(R,K);
A_MUA_global_all = cell(R,K);
M_MUA_k = cell(R,K);
A_MUA = cell(1,R);
A_MUA_final = cell(1,R);
MUA_bestrun_idx = zeros(1,R);
time_MUA = zeros(R,K);

fprintf('\n\nGMBUA')
for r = 1:R
    fprintf(['\nR = ' num2str(r) ':'])
    fprintf('\nK = ')
    for k = 1:K
        rng(k + (r-1)*K);
        fprintf([num2str(k) ', '])
        lambda   = 0.1;
        fraction = 0.3;
        beta     = 1;
        slic_size  = 7;
        slic_reg   = 0.001;
        % extract endmember bundles
        [groups{r,k}, bundle{r,k}] = batchvca2(X, P, bundle_nbr, percent);
        tic
        A_MUA_all{r,k} = genMUA_social(Xim,bundle{r,k},groups{r,k},A_init,slic_size,slic_reg,beta,lambda,...
            rho,maxiter_ADMM,'fractional',fraction,tol_a,0,vlfeat_path,0);
        time_MUA(r,k) = toc;
        [A_MUA_global_all{r,k}, S_MUA_global] = bundle2global(A_MUA_all{r,k},bundle{r,k},groups{r,k});
        M_MUA_k{r,k} = match_endmembers(A_MUA_global_all{r,k}, A_true_final);
        A_MUA_global_all{r,k} = A_MUA_global_all{r,k}(M_MUA_k{r,k}(:,1),:);
    end
    % find most representative abundance
    [A_MUA_final{r}, ~, MUA_bestrun_idx(r)] = find_most_consistent_run_abundances(A_MUA_global_all(r,:));    
    A_MUA{r} = A_MUA_all{r, MUA_bestrun_idx(r)};
end
fprintf('\n')

H_MUA = cell(1,R);
RMSE_MUA = zeros(1,R);
SAM_MUA = zeros(1,R);
SRE_MUA_k = zeros(R,K);
SRE_MUA = zeros(1,R);
for r=1:R
    H_MUA{1,r} = bundle{r, MUA_bestrun_idx(r)}*A_MUA{r}; 
    RMSE_px = sqrt(1/L*sum((H_MUA{r}-X).^2,1));
    RMSE_MUA(r) = mean(RMSE_px(:)); 
    SAM_px = zeros(N,1);
    for i = 1:N
        SAM_px(i) = 180/pi*real(acos((X(:,i)'*H_MUA{r}(:,i))...
            /(norm(X(:,i))*norm(H_MUA{r}(:,i)))));
    end
    SAM_MUA(r) = mean(SAM_px(:));

    for k=1:K
        SRE_MUA_k(r,k) = 20*log10(norm(A_true_final,'fro')/norm(A_MUA_global_all{r,k} - A_true_final,'fro'));
    end
    SRE_MUA(r) = SRE_MUA_k(r, MUA_bestrun_idx(r));
end
aux = sort(SRE_MUA);
SRE_MUA_median = aux(ceil(length(SRE_MUA)/2));
MUA_median_idx = find(SRE_MUA == SRE_MUA_median, 1, 'first');

%% Save results
disp('Saving to file...')
% save(strcat(path,image_name,'/',image_name,'_compare','_noise',num2str(noise),'_GMBUA_K',num2str(K),'_p',num2str(percent)),'-v7.3');

%%

figure, set(gcf,'color', 'white'),
plot(1:K_, SRE_FCLSU_k, '-o'), hold on,
plot(1:K_, SRE_collaborative_k, '-o'), hold on,
plot(1:K_, SRE_group_k, '-o'), hold on,
plot(1:K_, SRE_elitist_k, '-o'), hold on,
plot(1:K_, SRE_fractional_k, '-o'), hold on,
plot(1:R, SRE_MUA, '-o'), hold on,
scatter(FCLSU_median_idx, SRE_FCLSU_k(FCLSU_median_idx), 'red','*'),
scatter(FCLSU_bestrun_idx, SRE_FCLSU, 'blue','*'),
scatter(collaborative_median_idx, SRE_collaborative_k(collaborative_median_idx), 'red','*'),
scatter(collaborative_bestrun_idx, SRE_collaborative, 'blue','*'),
scatter(group_median_idx, SRE_group_k(group_median_idx), 'red','*'),
scatter(group_bestrun_idx, SRE_group, 'blue','*'),
scatter(elitist_median_idx, SRE_elitist_k(elitist_median_idx), 'red','*'),
scatter(elitist_bestrun_idx, SRE_elitist, 'blue','*'),
scatter(fractional_median_idx, SRE_fractional_k(fractional_median_idx), 'red','*'),
scatter(fractional_bestrun_idx, SRE_fractional, 'blue','*'),
scatter(MUA_median_idx, SRE_MUA(MUA_median_idx), 'red','*'),

xlabel('k'), xticks(1:K), xlim([0.5 K+0.5]), grid on,
ylabel('SRE (dB)'), ylim([0 15]),
legend({'FCLSU', 'collaborative', 'group', 'elitist', 'fractional', 'GMBUA', 'median', 'most consistent'},...
    'Location','northwest','NumColumns',2)
% exportgraphics(gcf,strcat(path,image_name,'/prints/',image_name,'_SRE_','K',num2str(K),'_',num2str(noise),'dB_p',num2str(percent),'_repl','.pdf'))

%% Display RMSEs, SAMs, SREs

fprintf('\n--- SAMs \n\n')
fprintf('FCLSU = %g \n', mean(SAM_FCLSU(:)))
fprintf('collaborative = %g \n', mean(SAM_collaborative(:)))
fprintf('group = %g \n', mean(SAM_group(:)))
fprintf('elitist = %g \n', mean(SAM_elitist(:)))
fprintf('fractional = %g \n', mean(SAM_fractional(:)))
fprintf('GMBUA = %g \n', SAM_MUA(r))

fprintf('\n--- RMSEs \n\n')
fprintf('FCLSU = %g \n', mean(RMSE_FCLSU(:)))
fprintf('collaborative = %g \n', mean(RMSE_collaborative(:)))
fprintf('group = %g \n', mean(RMSE_group(:)))
fprintf('elitist = %g \n', mean(RMSE_elitist(:)))
fprintf('fractional = %g \n', mean(RMSE_fractional(:)))
fprintf('GMBUA = %g \n', RMSE_MUA(r))

fprintf('\n--- SREs - Global Abundances \n\n')
fprintf('FCLSU = %g \n', SRE_FCLSU_median)
fprintf('collaborative = %g \n', SRE_collaborative_median)
fprintf('group = %g \n', SRE_group_median)
fprintf('elitist = %g \n', SRE_elitist_median)
fprintf('fractional = %g \n', SRE_fractional_median)
fprintf('GMBUA = %g \n', SRE_MUA(MUA_median_idx))

fprintf('\n--- Exec. - Time \n\n')
fprintf('FCLSU = %g \n', median(time_FCLSU))
fprintf('collaborative = %g \n', median(time_collaborative))
fprintf('group = %g \n', median(time_group))
fprintf('elitist = %g \n', median(time_elitist))
fprintf('fractional = %g \n', median(time_fractional))
fprintf('GMBUA = %g \n', median(time_MUA(1,:)))

%% Global abundance maps -- Most consistent

A_true_im = reshape(A_true_final',m,n,P);
A_FCLSU_im = reshape(A_FCLSU_final',m,n,P);
A_collaborative_im = reshape(A_collaborative_final',m,n,P);
A_group_im = reshape(A_group_final',m,n,P);
A_elitist_im = reshape(A_elitist_final',m,n,P);
A_fractional_im = reshape(A_fractional_final',m,n,P);
A_MUA_im = reshape(A_MUA_final{MUA_median_idx}',m,n,P);

vertical_spacing = 0.002;
horizontal_spacing = 0.002;
font_s = 13;
figure, set(gcf,'color', 'white')
for p = 1:P
    
    subaxis(7,P,p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_true_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('True','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_FCLSU_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('FCLSU','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,2*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_collaborative_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Collab.','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,3*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_group_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Group','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,4*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_elitist_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Elitist','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,5*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_fractional_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Fractional','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,6*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_MUA_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('GMBUA','fontname','times','fontsize',font_s)
        xlabel('Bundle 1','fontname','times','fontsize',font_s)
    elseif p == 2
        xlabel('Bundle 2','fontname','times','fontsize',font_s)
    elseif p == 3
        xlabel('Bundle 3','fontname','times','fontsize',font_s)
    elseif p == 4
        xlabel('Bundle 4','fontname','times','fontsize',font_s)
    elseif p == 5
        xlabel('Bundle 5','fontname','times','fontsize',font_s)
    elseif p == 6
        xlabel('Bundle 6','fontname','times','fontsize',font_s)
    elseif p == 7
        xlabel('Bundle 7','fontname','times','fontsize',font_s)
    elseif p == 8
        xlabel('Bundle 8','fontname','times','fontsize',font_s)
    elseif p == 9
        xlabel('Bundle 9','fontname','times','fontsize',font_s)
    end
end
% exportgraphics(gcf,strcat(path,image_name,'/prints/',image_name,'_compare_','K',num2str(K),'_',num2str(noise),'dB','.pdf'))

%% Global abundance maps -- Median SRE

A_true_im = reshape(A_true_final',m,n,P);
A_FCLSU_im = reshape(A_FCLSU_global_all{1,FCLSU_median_idx}',m,n,P);
A_collaborative_im =reshape(A_collaborative_global_all{1,collaborative_median_idx}',m,n,P);
A_group_im = reshape(A_group_global_all{1,group_median_idx}',m,n,P);
A_elitist_im = reshape(A_elitist_global_all{1,elitist_median_idx}',m,n,P);
A_fractional_im = reshape(A_fractional_global_all{1,fractional_median_idx}',m,n,P);
A_MUA_im = reshape(A_MUA_final{MUA_median_idx}',m,n,P);

vertical_spacing = 0.002;
horizontal_spacing = 0.002;
font_s = 13;
figure, set(gcf,'color', 'white')
for p = 1:P
    
    subaxis(7,P,p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_true_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('True','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_FCLSU_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('FCLSU','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,2*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_collaborative_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Collab.','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,3*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_group_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Group','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,4*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_elitist_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Elitist','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,5*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_fractional_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Fractional','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,6*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_MUA_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('GMBUA','fontname','times','fontsize',font_s)
        xlabel('Bundle 1','fontname','times','fontsize',font_s)
    elseif p == 2
        xlabel('Bundle 2','fontname','times','fontsize',font_s)
    elseif p == 3
        xlabel('Bundle 3','fontname','times','fontsize',font_s)
    elseif p == 4
        xlabel('Bundle 4','fontname','times','fontsize',font_s)
    elseif p == 5
        xlabel('Bundle 5','fontname','times','fontsize',font_s)
    elseif p == 6
        xlabel('Bundle 6','fontname','times','fontsize',font_s)
    elseif p == 7
        xlabel('Bundle 7','fontname','times','fontsize',font_s)
    elseif p == 8
        xlabel('Bundle 8','fontname','times','fontsize',font_s)
    elseif p == 9
        xlabel('Bundle 9','fontname','times','fontsize',font_s)
    end
end
% exportgraphics(gcf,strcat(path,image_name,'/prints/',image_name,'_compare_','K',num2str(K),'_',num2str(noise),'dB','.pdf'))


%% Reconstructed observed HI

X_true_im = reshape(H_true',m,n,L);
X_FCLSU_im = reshape(H_FCLSU',m,n,L);
X_collaborative_im =reshape(H_collaborative',m,n,L);
X_group_im = reshape(H_group',m,n,L);
X_elitist_im = reshape(H_elitist',m,n,L);
X_fractional_im = reshape(H_fractional',m,n,L);
X_MUA_im = reshape(H_MUA{MUA_median_idx}',m,n,L);
bands = [29 15 12];

vertical_spacing = 0.002;
horizontal_spacing = 0.002;
font_s = 13;
figure, set(gcf,'color', 'white')
subaxis(1,7,1,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
img = X_true_im(:,:,bands);
img = uint8(255*(img - min(img(:)))./(max(img(:))-min(img(:))));
imshow(img,[])
set(gca,'clim',[0,1])
xlabel('True','fontname','times','fontsize',font_s)

subaxis(1,7,2,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
img = X_FCLSU_im(:,:,bands);
img = uint8(255*(img - min(img(:)))./(max(img(:))-min(img(:))));
imshow(img,[])
set(gca,'clim',[0,1])
xlabel('FCLSU','fontname','times','fontsize',font_s)

subaxis(1,7,3,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
img = X_collaborative_im(:,:,bands);
img = uint8(255*(img - min(img(:)))./(max(img(:))-min(img(:))));
imshow(img,[])
set(gca,'clim',[0,1])
xlabel('Collab.','fontname','times','fontsize',font_s)

subaxis(1,7,4,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
img = X_group_im(:,:,bands);
img = uint8(255*(img - min(img(:)))./(max(img(:))-min(img(:))));
imshow(img,[])
set(gca,'clim',[0,1])
xlabel('Group','fontname','times','fontsize',font_s)

subaxis(1,7,5,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
img = X_elitist_im(:,:,bands);
img = uint8(255*(img - min(img(:)))./(max(img(:))-min(img(:))));
imshow(img,[])
set(gca,'clim',[0,1])
xlabel('Elitist','fontname','times','fontsize',font_s)

subaxis(1,7,6,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
img = X_fractional_im(:,:,bands);
img = uint8(255*(img - min(img(:)))./(max(img(:))-min(img(:))));
imshow(img,[])
set(gca,'clim',[0,1])
xlabel('Fractional','fontname','times','fontsize',font_s)

subaxis(1,7,7,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
img = X_MUA_im(:,:,bands);
img = uint8(255*(img - min(img(:)))./(max(img(:))-min(img(:))));
imshow(img,[])
xlabel('GMBUA','fontname','times','fontsize',font_s)
set(gca,'clim',[0,1])
% export_fig(strcat(path,image_name,'/prints/',image_name,'_',test_name,'_compare_','frac_',num2str(noise),'dB'),'-pdf','-pdf');

%% Plot bundles abundances per pixel

vertical_spacing = 0.002;
horizontal_spacing = 0.002;
font_s = 11;
figure, set(gcf,'color', 'white')

subaxis(1,7,1,'spacinghorizontal',horizontal_spacing,'spacingvertical',vertical_spacing);
imagesc(A_true_final),
set(gca,'fontname','times','fontsize',font_s)
xticks([]), yticks([]), axis square
xlabel('True','fontname','times','fontsize',font_s)
colormap jet

subaxis(1,7,2,'spacinghorizontal',horizontal_spacing,'spacingvertical',vertical_spacing);
imagesc(A_FCLSU_final)
set(gca,'fontname','times','fontsize',font_s)
xticks([]), yticks([]), axis square
xlabel('FCLSU','fontname','times','fontsize',font_s)
colormap jet

subaxis(1,7,3,'spacinghorizontal',horizontal_spacing,'spacingvertical',vertical_spacing);
imagesc(A_collaborative_final)
set(gca,'fontname','times','fontsize',font_s)
xticks([]), yticks([]), axis square
xlabel('Collab.','fontname','times','fontsize',font_s)
colormap jet

subaxis(1,7,4,'spacinghorizontal',horizontal_spacing,'spacingvertical',vertical_spacing);
imagesc(A_group_final)
set(gca,'fontname','times','fontsize',font_s)
xticks([]), yticks([]), axis square
xlabel('Group','fontname','times','fontsize',font_s)
colormap jet

subaxis(1,7,5,'spacinghorizontal',horizontal_spacing,'spacingvertical',vertical_spacing);
imagesc(A_elitist_final)
set(gca,'fontname','times','fontsize',font_s)
xticks([]), yticks([]), axis square
xlabel('Elitist','fontname','times','fontsize',font_s)
colormap jet

subaxis(1,7,6,'spacinghorizontal',horizontal_spacing,'spacingvertical',vertical_spacing);
imagesc(A_fractional_final)
set(gca,'fontname','times','fontsize',font_s)
xticks([]), yticks([]), axis square
xlabel('Fractional','fontname','times','fontsize',font_s)
colormap jet

subaxis(1,7,7,'spacinghorizontal',horizontal_spacing,'spacingvertical',vertical_spacing);
imagesc(A_MUA_final{MUA_median_idx})
set(gca,'fontname','times','fontsize',font_s)
xticks([]), yticks([]), axis square
xlabel('GMBUA','fontname','times','fontsize',font_s)
colormap jet

%% Plot endmember bundles

% Ground-truth
plot_bundles(bundle_true, groups_true);
% Extracted in most consistent run
plot_bundles(bundle{MUA_median_idx,MUA_bestrun_idx(MUA_median_idx)}, groups{MUA_median_idx,MUA_bestrun_idx(MUA_median_idx)}, M_MUA_k{MUA_median_idx,MUA_bestrun_idx(MUA_median_idx)});
% exportgraphics(gcf,strcat(path,image_name,'/prints/',image_name,'_bundles_','K',num2str(K),'_',num2str(noise),'dB','.pdf'))
