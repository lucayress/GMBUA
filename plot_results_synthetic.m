%% Show Results
clc, clearvars, close all
tic

path = 'test_data/';
image_name = 'synthetic4';

load(strcat(path,image_name,'/results/','synthetic4_compare_all_noise10_R30'));
% load(strcat(path,image_name,'/results/','synthetic4_compare_all_noise20_R30'));
% load(strcat(path,image_name,'/results/','synthetic4_compare_all_noise30_R30'));

%% Display RMSEs, SAMs, SREs

fprintf('\n--- SAMs \n\n')
fprintf('FCLSU = %g \n', mean(SAM_FCLSU(:)))
fprintf('collaborative = %g \n', mean(SAM_collaborative(:)))
fprintf('group = %g \n', mean(SAM_group(:)))
fprintf('elitist = %g \n', mean(SAM_elitist(:)))
fprintf('fractional = %g \n', mean(SAM_fractional(:)))
fprintf('SUnCNN = %g \n', mean(SAM_SUnCNN(:)))
fprintf('GMBUA = %g \n', SAM_MUA(r))

fprintf('\n--- RMSEs \n\n')
fprintf('FCLSU = %g \n', mean(RMSE_FCLSU(:)))
fprintf('collaborative = %g \n', mean(RMSE_collaborative(:)))
fprintf('group = %g \n', mean(RMSE_group(:)))
fprintf('elitist = %g \n', mean(RMSE_elitist(:)))
fprintf('fractional = %g \n', mean(RMSE_fractional(:)))
fprintf('SUnCNN = %g \n', mean(RMSE_SUnCNN(:)))
fprintf('GMBUA = %g \n', RMSE_MUA(r))

fprintf('\n--- SREs - Global Abundances \n\n')
fprintf('FCLSU = %g \n', SRE_FCLSU_median)
fprintf('collaborative = %g \n', SRE_collaborative_median)
fprintf('group = %g \n', SRE_group_median)
fprintf('elitist = %g \n', SRE_elitist_median)
fprintf('fractional = %g \n', SRE_fractional_median)
fprintf('SUnCNN = %g \n', SRE_SUnCNN_median)
fprintf('GMBUA = %g \n', SRE_MUA(MUA_median_idx))

SRE_FCLSU_X = 20*log10(norm(H_true,'fro')/norm(H_FCLSU - H_true,'fro'));
SRE_collaborative_X =  20*log10(norm(H_true,'fro')/norm(H_collaborative - H_true,'fro'));
SRE_group_X = 20*log10(norm(H_true,'fro')/norm(H_group - H_true,'fro'));
SRE_elitist_X = 20*log10(norm(H_true,'fro')/norm(H_elitist - H_true,'fro'));
SRE_fractional_X = 20*log10(norm(H_true,'fro')/norm(H_fractional - H_true,'fro'));
SRE_MUA_X = 20*log10(norm(H_true,'fro')/norm(H_MUA{1,MUA_median_idx} - H_true,'fro'));

fprintf('\n--- SREs - Observed Image \n\n')
fprintf('FCLSU = %g \n', SRE_FCLSU_X)
fprintf('collaborative = %g \n', SRE_collaborative_X)
fprintf('group = %g \n', SRE_group_X)
fprintf('elitist = %g \n', SRE_elitist_X)
fprintf('fractional = %g \n', SRE_fractional_X)
fprintf('SUnCNN = %g \n', SRE_SUnCNN_X)
fprintf('GMBUA = %g \n', SRE_MUA_X)

%% Global abundance maps -- Most consistent

A_true_im = reshape(A_true_final',m,n,P);
A_FCLSU_im = reshape(A_FCLSU_final',m,n,P);
A_collaborative_im = reshape(A_collaborative_final',m,n,P);
A_group_im = reshape(A_group_final',m,n,P);
A_elitist_im = reshape(A_elitist_final',m,n,P);
A_fractional_im = reshape(A_fractional_final',m,n,P);
A_SUnCNN_im = reshape(A_SUnCNN_final',m,n,P);
A_MUA_im = reshape(A_MUA_final{MUA_median_idx}',m,n,P);

vertical_spacing = 0.002;
horizontal_spacing = 0.002;
font_s = 13;
figure, set(gcf,'color', 'white')
for p = 1:P
    
    subaxis(8,P,p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_true_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('GT','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(8,P,P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_FCLSU_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('FCLSU','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(8,P,2*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_collaborative_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Collab.','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(8,P,3*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_group_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Group','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(8,P,4*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_elitist_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Elitist','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(8,P,5*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_fractional_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Fractional','fontname','times','fontsize',font_s)
    end
    colormap jet

    subaxis(8,P,6*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_SUnCNN_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('SUnCNN','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(8,P,7*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_MUA_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('GMBUA','fontname','times','fontsize',font_s)
        xlabel('Endm. 1','fontname','times','fontsize',font_s)
    elseif p == 2
        xlabel('Endm. 2','fontname','times','fontsize',font_s)
    elseif p == 3
        xlabel('Endm. 3','fontname','times','fontsize',font_s)
    elseif p == 4
        xlabel('Endm. 4','fontname','times','fontsize',font_s)
    elseif p == 5
        xlabel('Endm. 5','fontname','times','fontsize',font_s)
    elseif p == 6
        xlabel('Endm. 6','fontname','times','fontsize',font_s)
    elseif p == 7
        xlabel('Endm. 7','fontname','times','fontsize',font_s)
    elseif p == 8
        xlabel('Endm. 8','fontname','times','fontsize',font_s)
    elseif p == 9
        xlabel('Endm. 9','fontname','times','fontsize',font_s)
    end
end
% Create colorbar
colorbar('Position',...
    [0.893337956919058 0.100701943844492 0.00833333333333341 0.112041036717063]);

% exportgraphics(gcf,strcat(path,image_name,'/prints/',image_name,'_compare_','K',num2str(K),'_',num2str(SNR),'dB','.pdf'))

%% Global abundance maps -- Median SRE

A_true_im = reshape(A_true_final',m,n,P);
A_FCLSU_im = reshape(A_FCLSU_global_all{1,FCLSU_median_idx}',m,n,P);
A_collaborative_im =reshape(A_collaborative_global_all{1,collaborative_median_idx}',m,n,P);
A_group_im = reshape(A_group_global_all{1,group_median_idx}',m,n,P);
A_elitist_im = reshape(A_elitist_global_all{1,elitist_median_idx}',m,n,P);
A_fractional_im = reshape(A_fractional_global_all{1,fractional_median_idx}',m,n,P);
A_SUnCNN_im = reshape(A_SUnCNN_global_all{1,SUnCNN_median_idx}',m,n,P);
A_MUA_im = reshape(A_MUA_final{MUA_median_idx}',m,n,P);

vertical_spacing = 0.002;
horizontal_spacing = 0.002;
font_s = 13;
figure, set(gcf,'color', 'white')
for p = 1:P
    
    subaxis(8,P,p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_true_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('GT','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(8,P,P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_FCLSU_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('FCLSU','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(8,P,2*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_collaborative_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Collab.','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(8,P,3*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_group_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Group','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(8,P,4*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_elitist_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Elitist','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(8,P,5*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_fractional_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Fractional','fontname','times','fontsize',font_s)
    end
    colormap jet

    subaxis(8,P,6*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_SUnCNN_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('SUnCNN','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(8,P,7*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_MUA_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('GMBUA','fontname','times','fontsize',font_s)
        xlabel('Endm. 1','fontname','times','fontsize',font_s)
    elseif p == 2
        xlabel('Endm. 2','fontname','times','fontsize',font_s)
    elseif p == 3
        xlabel('Endm. 3','fontname','times','fontsize',font_s)
    elseif p == 4
        xlabel('Endm. 4','fontname','times','fontsize',font_s)
    elseif p == 5
        xlabel('Endm. 5','fontname','times','fontsize',font_s)
    elseif p == 6
        xlabel('Endm. 6','fontname','times','fontsize',font_s)
    elseif p == 7
        xlabel('Endm. 7','fontname','times','fontsize',font_s)
    elseif p == 8
        xlabel('Endm. 8','fontname','times','fontsize',font_s)
    elseif p == 9
        xlabel('Endm. 9','fontname','times','fontsize',font_s)
    end
end
% Create colorbar
colorbar('Position',...
    [0.893337956919058 0.100701943844492 0.00833333333333341 0.112041036717063]);
% exportgraphics(gcf,strcat(path,image_name,'/prints/','synthetic4_compare_all_noise30_R30','.pdf'))

%% Global abundance maps -- Median SRE - 5 algorithms

vertical_spacing = 0.002;
horizontal_spacing = 0.002;
font_s = 10;
figure, set(gcf,'color', 'white')
for p = 1:P
    
    subaxis(5,P,p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_true_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('GT','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(5,P,P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_elitist_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Elitist','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(5,P,2*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_fractional_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Fractional','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(5,P,3*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_SUnCNN_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('SUnCNN','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(5,P,4*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_MUA_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('GMBUA','fontname','times','fontsize',font_s)
        xlabel('Endm. 1','fontname','times','fontsize',font_s)
    elseif p == 2
        xlabel('Endm. 2','fontname','times','fontsize',font_s)
    elseif p == 3
        xlabel('Endm. 3','fontname','times','fontsize',font_s)
    elseif p == 4
        xlabel('Endm. 4','fontname','times','fontsize',font_s)
    elseif p == 5
        xlabel('Endm. 5','fontname','times','fontsize',font_s)
    elseif p == 6
        xlabel('Endm. 6','fontname','times','fontsize',font_s)
    elseif p == 7
        xlabel('Endm. 7','fontname','times','fontsize',font_s)
    elseif p == 8
        xlabel('Endm. 8','fontname','times','fontsize',font_s)
    elseif p == 9
        xlabel('Endm. 9','fontname','times','fontsize',font_s)
    end
end
% Create colorbar
colorbar('Position',...
    [0.900835073068894 0.100709606986899 0.0137743207857127 0.156932314410481]);
% exportgraphics(gcf,strcat(path,image_name,'/prints/','synthetic4_compare_all_noise30_R30_5algs','.pdf'))

%% Plot reconstructed observed HI

X_true_im = reshape(H_true',m,n,L);
X_FCLSU_im = reshape(H_FCLSU',m,n,L);
X_collaborative_im =reshape(H_collaborative',m,n,L);
X_group_im = reshape(H_group',m,n,L);
X_elitist_im = reshape(H_elitist',m,n,L);
X_fractional_im = reshape(H_fractional',m,n,L);
X_SUnCNN_im = reshape(H_SUnCNN',m,n,L);
X_MUA_im = reshape(H_MUA{MUA_median_idx}',m,n,L);
bands = [29 15 12];

vertical_spacing = 0.002;
horizontal_spacing = 0.002;
font_s = 13;
figure, set(gcf,'color', 'white')
subaxis(1,8,1,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
img = X_true_im(:,:,bands);
img = uint8(255*(img - min(img(:)))./(max(img(:))-min(img(:))));
imshow(img,[])
set(gca,'clim',[0,1])
xlabel('GT','fontname','times','fontsize',font_s)

subaxis(1,8,2,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
img = X_FCLSU_im(:,:,bands);
img = uint8(255*(img - min(img(:)))./(max(img(:))-min(img(:))));
imshow(img,[])
set(gca,'clim',[0,1])
xlabel('FCLSU','fontname','times','fontsize',font_s)

subaxis(1,8,3,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
img = X_collaborative_im(:,:,bands);
img = uint8(255*(img - min(img(:)))./(max(img(:))-min(img(:))));
imshow(img,[])
set(gca,'clim',[0,1])
xlabel('Collab.','fontname','times','fontsize',font_s)

subaxis(1,8,4,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
img = X_group_im(:,:,bands);
img = uint8(255*(img - min(img(:)))./(max(img(:))-min(img(:))));
imshow(img,[])
set(gca,'clim',[0,1])
xlabel('Group','fontname','times','fontsize',font_s)

subaxis(1,8,5,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
img = X_elitist_im(:,:,bands);
img = uint8(255*(img - min(img(:)))./(max(img(:))-min(img(:))));
imshow(img,[])
set(gca,'clim',[0,1])
xlabel('Elitist','fontname','times','fontsize',font_s)

subaxis(1,8,6,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
img = X_fractional_im(:,:,bands);
img = uint8(255*(img - min(img(:)))./(max(img(:))-min(img(:))));
imshow(img,[])
set(gca,'clim',[0,1])
xlabel('Fractional','fontname','times','fontsize',font_s)

subaxis(1,8,7,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
img = X_SUnCNN_im(:,:,bands);
img = uint8(255*(img - min(img(:)))./(max(img(:))-min(img(:))));
imshow(img,[])
set(gca,'clim',[0,1])
xlabel('SUnCNN','fontname','times','fontsize',font_s)

subaxis(1,8,8,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
img = X_MUA_im(:,:,bands);
img = uint8(255*(img - min(img(:)))./(max(img(:))-min(img(:))));
imshow(img,[])
xlabel('GMBUA','fontname','times','fontsize',font_s)
set(gca,'clim',[0,1])
% export_fig(strcat(path,image_name,'/prints/',image_name,'_',test_name,'_compare_','frac_',num2str(SNR),'dB'),'-pdf','-pdf');

%% Plot bundles abundances per pixel

vertical_spacing = 0.002;
horizontal_spacing = 0.002;
font_s = 11;
figure, set(gcf,'color', 'white')

subaxis(1,8,1,'spacinghorizontal',horizontal_spacing,'spacingvertical',vertical_spacing);
imagesc(A_true_final),
set(gca,'fontname','times','fontsize',font_s)
xticks([]), yticks([]), axis square
xlabel('GT','fontname','times','fontsize',font_s)
colormap jet

subaxis(1,8,2,'spacinghorizontal',horizontal_spacing,'spacingvertical',vertical_spacing);
imagesc(A_FCLSU_final)
set(gca,'fontname','times','fontsize',font_s)
xticks([]), yticks([]), axis square
xlabel('FCLSU','fontname','times','fontsize',font_s)
colormap jet

subaxis(1,8,3,'spacinghorizontal',horizontal_spacing,'spacingvertical',vertical_spacing);
imagesc(A_collaborative_final)
set(gca,'fontname','times','fontsize',font_s)
xticks([]), yticks([]), axis square
xlabel('Collab.','fontname','times','fontsize',font_s)
colormap jet

subaxis(1,8,4,'spacinghorizontal',horizontal_spacing,'spacingvertical',vertical_spacing);
imagesc(A_group_final)
set(gca,'fontname','times','fontsize',font_s)
xticks([]), yticks([]), axis square
xlabel('Group','fontname','times','fontsize',font_s)
colormap jet

subaxis(1,8,5,'spacinghorizontal',horizontal_spacing,'spacingvertical',vertical_spacing);
imagesc(A_elitist_final)
set(gca,'fontname','times','fontsize',font_s)
xticks([]), yticks([]), axis square
xlabel('Elitist','fontname','times','fontsize',font_s)
colormap jet

subaxis(1,8,6,'spacinghorizontal',horizontal_spacing,'spacingvertical',vertical_spacing);
imagesc(A_fractional_final)
set(gca,'fontname','times','fontsize',font_s)
xticks([]), yticks([]), axis square
xlabel('Fractional','fontname','times','fontsize',font_s)
colormap jet

subaxis(1,8,7,'spacinghorizontal',horizontal_spacing,'spacingvertical',vertical_spacing);
imagesc(A_SUnCNN_final)
set(gca,'fontname','times','fontsize',font_s)
xticks([]), yticks([]), axis square
xlabel('SUnCNN','fontname','times','fontsize',font_s)
colormap jet

subaxis(1,8,8,'spacinghorizontal',horizontal_spacing,'spacingvertical',vertical_spacing);
imagesc(A_MUA_final{MUA_median_idx})
set(gca,'fontname','times','fontsize',font_s)
xticks([]), yticks([]), axis square
xlabel('GMBUA','fontname','times','fontsize',font_s)
colormap jet

%% SRE plot

figure, set(gcf,'color', 'white'),
plot(1:K, SRE_FCLSU_k, '-o'), hold on,
plot(1:K, SRE_collaborative_k, '-o'), hold on,
plot(1:K, SRE_group_k, '-o'), hold on,
plot(1:K, SRE_elitist_k, '-o'), hold on,
plot(1:K, SRE_fractional_k, '-o'), hold on,
plot(1:K, SRE_SUnCNN_k, '-o'), hold on,
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
scatter(SUnCNN_median_idx, SRE_SUnCNN_k(SUnCNN_median_idx), 'red','*'),
scatter(SUnCNN_bestrun_idx, SRE_SUnCNN, 'blue','*'),
scatter(MUA_median_idx, SRE_MUA(MUA_median_idx), 'red','*'),

xlabel('k'), xticks(1:K), xlim([0.5 K+0.5]), grid on,
ylabel('SRE (dB)'), ylim([0 15]),
legend({'FCLSU', 'collaborative', 'group', 'elitist', 'fractional', 'SUnCNN', 'GMBUA', 'median', 'most consistent'},...
    'Location','northwest','NumColumns',2)
% exportgraphics(gcf,strcat(path,image_name,'/prints/',image_name,'_SRE_','K',num2str(K),'_',num2str(SNR),'dB_p',num2str(percent),'_repl','.pdf'))

%% Plot endmember bundles

% Ground-truth
plot_bundles(bundle_true, groups_true);
% Extracted in most consistent run
plot_bundles(bundle{MUA_median_idx,MUA_bestrun_idx(MUA_median_idx)}, groups{MUA_median_idx,MUA_bestrun_idx(MUA_median_idx)}, M_MUA_k{MUA_median_idx,MUA_bestrun_idx(MUA_median_idx)});
% exportgraphics(gcf,strcat(path,image_name,'/prints/',image_name,'_bundles_','K',num2str(K),'_',num2str(SNR),'dB','.pdf'))
