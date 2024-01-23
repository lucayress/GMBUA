%% Show Results
clc, clearvars, close all
tic

path = 'real_data/';
image_name = 'cuprite';

load(strcat(path,image_name,'/results/','cuprite_compare_all_noise0_R30.mat'));

%% Display RMSEs, SAMs, SREs

fprintf('\n--- SAMs \n\n')
fprintf('FCLSU = %g \n', mean(SAM_FCLSU(:)))
fprintf('collaborative = %g \n', mean(SAM_collaborative(:)))
fprintf('group = %g \n', mean(SAM_group(:)))
fprintf('elitist = %g \n', mean(SAM_elitist(:)))
fprintf('fractional = %g \n', mean(SAM_fractional(:)))
fprintf('SUnCNN = %g \n', mean(SAM_SUnCNN(:)))
fprintf('GMBUA = %g \n', SAM_MUA(MUA_median_idx))

fprintf('\n--- RMSEs \n\n')
fprintf('FCLSU = %g \n', mean(RMSE_FCLSU(:)))
fprintf('collaborative = %g \n', mean(RMSE_collaborative(:)))
fprintf('group = %g \n', mean(RMSE_group(:)))
fprintf('elitist = %g \n', mean(RMSE_elitist(:)))
fprintf('fractional = %g \n', mean(RMSE_fractional(:)))
fprintf('SUnCNN = %g \n', mean(RMSE_SUnCNN(:)))
fprintf('GMBUA = %g \n', RMSE_MUA(MUA_median_idx))

fprintf('\n--- SREs - Observed Image \n\n')
fprintf('FCLSU = %g \n', SRE_FCLSU_median)
fprintf('collaborative = %g \n', SRE_collaborative_median)
fprintf('group = %g \n', SRE_group_median)
fprintf('elitist = %g \n', SRE_elitist_median)
fprintf('fractional = %g \n', SRE_fractional_median)
fprintf('SUnCNN = %g \n', SRE_SUnCNN_median)
fprintf('GMBUA = %g \n', SRE_MUA_median)

fprintf('\n--- Exec. - Time \n\n')
fprintf('SUnCNN = %g \n', median(time_SUnCNN))

%% Line plot SRE(Y)
figure, set(gcf,'color', 'white'),
plot(1:K, SRE_FCLSU_k, '-o'), hold on,
plot(1:K, SRE_collaborative_k, '-o'), hold on,
plot(1:K, SRE_group_k, '-o'), hold on,
plot(1:K, SRE_elitist_k, '-o'), hold on,
plot(1:K, SRE_fractional_k, '-o'), hold on,
plot(1:K, SRE_SUnCNN_k, '-o'), hold on,
plot(1:length(seed_R), SRE_MUA, '-o'), hold on,
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
ylabel('SRE (dB)'), ylim([30 45]),
legend({'FCLSU', 'collaborative', 'group', 'elitist', 'fractional', 'SUnCNN', 'GMBUA', 'median', 'most consistent'},...
    'Location','northwest','NumColumns',2)
% exportgraphics(gcf,strcat(path,image_name,'/prints/',image_name,'_compare_noise',num2str(noise),'_R',num2str(length(seed_R)),'_K',num2str(K),'.pdf'))

%% Global abundance maps -- Most consistent

A_FCLSU_im = reshape(A_FCLSU_global_all{1,FCLSU_bestrun_idx}',m,n,P);
A_collaborative_im =reshape(A_collaborative_global_all{1,collaborative_bestrun_idx}',m,n,P);
A_group_im = reshape(A_group_global_all{1,group_bestrun_idx}',m,n,P);
A_elitist_im = reshape(A_elitist_global_all{1,elitist_bestrun_idx}',m,n,P);
A_fractional_im = reshape(A_fractional_global_all{1,fractional_bestrun_idx}',m,n,P);
A_SUnCNN_im = reshape(A_SUnCNN_global_all{1,SUnCNN_bestrun_idx}',m,n,P);
A_MUA_im = reshape(A_MUA_final{MUA_median_idx}',m,n,P);

vertical_spacing = 0.002;
horizontal_spacing = 0.002;
font_s = 13;
figure, set(gcf,'color', 'white')
for p = 1:P
    subaxis(7,P,p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_FCLSU_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('FCLSU','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_collaborative_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Collab.','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,2*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_group_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Group','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,3*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_elitist_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Elitist','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,4*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_fractional_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Fractional','fontname','times','fontsize',font_s)
    end
    colormap jet

    subaxis(7,P,5*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_SUnCNN_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('SUnCNN','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,6*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_MUA_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('GMBUA','fontname','times','fontsize',font_s)
        xlabel('Em 1','fontname','times','fontsize',font_s)
    elseif p == 2
        xlabel('Em 2','fontname','times','fontsize',font_s)
    elseif p == 3
        xlabel('Em 3','fontname','times','fontsize',font_s)
    elseif p == 4
        xlabel('Em 4','fontname','times','fontsize',font_s)
    elseif p == 5
        xlabel('Em 5','fontname','times','fontsize',font_s)
    elseif p == 6
        xlabel('Em 6','fontname','times','fontsize',font_s)
    elseif p == 7
        xlabel('Em 7','fontname','times','fontsize',font_s)
    elseif p == 8
        xlabel('Em 8','fontname','times','fontsize',font_s)
    elseif p == 9
        xlabel('Em 9','fontname','times','fontsize',font_s)
    elseif p == 10
        xlabel('Em 10','fontname','times','fontsize',font_s)
    elseif p == 11
        xlabel('Em 11','fontname','times','fontsize',font_s)
    elseif p == 12
        xlabel('Em 12','fontname','times','fontsize',font_s)
    elseif p == 13
        xlabel('Em 13','fontname','times','fontsize',font_s)
    elseif p == 14
        xlabel('Em 14','fontname','times','fontsize',font_s)
    end
end

%% Global abundance maps -- Median SRE

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
    subaxis(7,P,p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_FCLSU_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('FCLSU','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_collaborative_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Collab.','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,2*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_group_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Group','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,3*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_elitist_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Elitist','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,4*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_fractional_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Fractional','fontname','times','fontsize',font_s)
    end
    colormap jet

    subaxis(7,P,5*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_SUnCNN_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('SUnCNN','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,P,6*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_MUA_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('GMBUA','fontname','times','fontsize',font_s)
        xlabel('Em 1','fontname','times','fontsize',font_s)
    elseif p == 2
        xlabel('Em 2','fontname','times','fontsize',font_s)
    elseif p == 3
        xlabel('Em 3','fontname','times','fontsize',font_s)
    elseif p == 4
        xlabel('Em 4','fontname','times','fontsize',font_s)
    elseif p == 5
        xlabel('Em 5','fontname','times','fontsize',font_s)
    elseif p == 6
        xlabel('Em 6','fontname','times','fontsize',font_s)
    elseif p == 7
        xlabel('Em 7','fontname','times','fontsize',font_s)
    elseif p == 8
        xlabel('Em 8','fontname','times','fontsize',font_s)
    elseif p == 9
        xlabel('Em 9','fontname','times','fontsize',font_s)
    elseif p == 10
        xlabel('Em 10','fontname','times','fontsize',font_s)
    elseif p == 11
        xlabel('Em 11','fontname','times','fontsize',font_s)
    elseif p == 12
        xlabel('Em 12','fontname','times','fontsize',font_s)
    elseif p == 13
        xlabel('Em 13','fontname','times','fontsize',font_s)
    elseif p == 14
        xlabel('Em 14','fontname','times','fontsize',font_s)
    end
end

%% Global abundance maps -- Most consistent -- select

ems = [7 10 6 9 1 14];

A_FCLSU_im = reshape(A_FCLSU_global_all{1,FCLSU_bestrun_idx}',m,n,P);
A_collaborative_im =reshape(A_collaborative_global_all{1,collaborative_bestrun_idx}',m,n,P);
A_group_im = reshape(A_group_global_all{1,group_bestrun_idx}',m,n,P);
A_elitist_im = reshape(A_elitist_global_all{1,elitist_bestrun_idx}',m,n,P);
A_fractional_im = reshape(A_fractional_global_all{1,fractional_bestrun_idx}',m,n,P);
A_SUnCNN_im = reshape(A_SUnCNN_global_all{1,SUnCNN_bestrun_idx}',m,n,P);
A_MUA_im = reshape(A_MUA_final{MUA_median_idx}',m,n,P);

vertical_spacing = 0.002;
horizontal_spacing = 0.002;
font_s = 13;
figure, set(gcf,'color', 'white')
for p = 1:length(ems)
    subaxis(6,length(ems),p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_FCLSU_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('FCLSU','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,length(ems),length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_collaborative_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Collab.','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,length(ems),2*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_group_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Group','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,length(ems),3*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_elitist_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Elitist','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,length(ems),4*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_fractional_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Fractional','fontname','times','fontsize',font_s)
    end
    colormap jet

    subaxis(7,length(ems),5*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_SUnCNN_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('SUnCNN','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,length(ems),6*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_MUA_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('GMBUA','fontname','times','fontsize',font_s)
        xlabel('Alunite','fontname','times','fontsize',font_s)
    elseif p == 2
        xlabel('Muscovite','fontname','times','fontsize',font_s)
    elseif p == 3
        xlabel('Kaolinite','fontname','times','fontsize',font_s)
    elseif p == 4
        xlabel('Sphene','fontname','times','fontsize',font_s)
    elseif p == 5
        xlabel('Buddingtonite','fontname','times','fontsize',font_s)
    elseif p == 6
        xlabel('Chalcedony','fontname','times','fontsize',font_s)
    end
end

%% Global abundance maps -- Median SRE -- select

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
for p = 1:length(ems)
    subaxis(7,length(ems),p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_FCLSU_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('FCLSU','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,length(ems),length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_collaborative_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Collab.','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,length(ems),2*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_group_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Group','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,length(ems),3*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_elitist_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Elitist','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,length(ems),4*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_fractional_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Fractional','fontname','times','fontsize',font_s)
    end
    colormap jet

    subaxis(7,length(ems),5*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_SUnCNN_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('SUnCNN','fontname','times','fontsize',font_s)
    end
    colormap jet
    
    subaxis(7,length(ems),6*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_MUA_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('GMBUA','fontname','times','fontsize',font_s)
        xlabel('Alunite','fontname','times','fontsize',font_s)
    elseif p == 2
        xlabel('Muscovite','fontname','times','fontsize',font_s)
    elseif p == 3
        xlabel('Kaolinite','fontname','times','fontsize',font_s)
    elseif p == 4
        xlabel('Sphene','fontname','times','fontsize',font_s)
    elseif p == 5
        xlabel('Buddingtonite','fontname','times','fontsize',font_s)
    elseif p == 6
        xlabel('Chalcedony','fontname','times','fontsize',font_s)
    end
end
% exportgraphics(gcf,strcat(path,image_name,'/prints/',image_name,'_compare_','_R',num2str(length(seed_R)),'_K',num2str(K),'_medianSRE','.pdf'))


%% Global abundance maps -- Manual select - GMBUA vs. other algorithms random pick

rng(40)
k1 = randsample(K,1);
k2 = randsample(K,1);
k3 = randsample(K,1);
k4 = randsample(K,1);
k5 = randsample(K,1);
k6 = randsample(K,1);
k7 = randsample(K,1);

A_FCLSU_im = reshape(A_FCLSU_global_all{1,k1}',m,n,P);
A_collaborative_im =reshape(A_collaborative_global_all{1,k2}',m,n,P);
A_group_im = reshape(A_group_global_all{1,k3}',m,n,P);
A_elitist_im = reshape(A_elitist_global_all{1,k4}',m,n,P);
A_fractional_im = reshape(A_fractional_global_all{1,k5}',m,n,P);
A_SUnCNN_im = reshape(A_SUnCNN_global_all{1,k7}',m,n,P);
A_MUA_im = reshape(A_MUA_global_all{r,k6}',m,n,P);

vertical_spacing = 0.002;
horizontal_spacing = 0.002;
font_s = 13;
figure, set(gcf,'color', 'white')
for p = 1:length(ems)
    subaxis(7,length(ems),p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_FCLSU_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('FCLSU','fontname','times','fontsize',font_s)
    end
    colormap jet

    subaxis(7,length(ems),length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_collaborative_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Collab.','fontname','times','fontsize',font_s)
    end
    colormap jet

    subaxis(7,length(ems),2*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_group_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Group','fontname','times','fontsize',font_s)
    end
    colormap jet

    subaxis(7,length(ems),3*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_elitist_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Elitist','fontname','times','fontsize',font_s)
    end
    colormap jet

    subaxis(7,length(ems),4*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_fractional_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Fractional','fontname','times','fontsize',font_s)
    end
    colormap jet

    subaxis(7,length(ems),5*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_SUnCNN_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('SUnCNN','fontname','times','fontsize',font_s)
    end
    colormap jet

    subaxis(7,length(ems),6*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_MUA_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square, axis square
    if p == 1
        ylabel('GMBUA','fontname','times','fontsize',font_s)
        xlabel('Alunite','fontname','times','fontsize',font_s)
    elseif p == 2
        xlabel('Muscovite','fontname','times','fontsize',font_s)
    elseif p == 3
        xlabel('Kaolinite','fontname','times','fontsize',font_s)
    elseif p == 4
        xlabel('Sphene','fontname','times','fontsize',font_s)
    elseif p == 5
        xlabel('Buddingtonite','fontname','times','fontsize',font_s)
    elseif p == 6
        xlabel('Chalcedony','fontname','times','fontsize',font_s)
    end
end
% exportgraphics(gcf,strcat(path,image_name,'/prints/','cuprite_compare_all_noise0_R30','.pdf'))

% Only 4 algorithms
figure, set(gcf,'color', 'white')
for p = 1:length(ems)
    subaxis(4,length(ems),p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_elitist_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Elitist','fontname','times','fontsize',font_s)
    end
    colormap jet
        
    subaxis(4,length(ems),length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_fractional_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('Fractional','fontname','times','fontsize',font_s)
    end
    colormap jet

    subaxis(4,length(ems),2*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_SUnCNN_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square
    if p == 1
        ylabel('SUnCNN','fontname','times','fontsize',font_s)
    end
    colormap jet

    subaxis(4,length(ems),3*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
    imshow(A_MUA_im(:,:,ems(p)),[],'colormap', jet)
    set(gca,'clim',[0,1]), axis square, axis square
    if p == 1
        ylabel('GMBUA','fontname','times','fontsize',font_s)
        xlabel('Alunite','fontname','times','fontsize',font_s)
    elseif p == 2
        xlabel('Muscovite','fontname','times','fontsize',font_s)
    elseif p == 3
        xlabel('Kaolinite','fontname','times','fontsize',font_s)
    elseif p == 4
        xlabel('Sphene','fontname','times','fontsize',font_s)
    elseif p == 5
        xlabel('Buddingtonite','fontname','times','fontsize',font_s)
    elseif p == 6
        xlabel('Chalcedony','fontname','times','fontsize',font_s)
    end
end
% Create colorbar
colorbar('Position',...
    [0.903880719819764 0.105633802816901 0.00720473052665682 0.188898825049061]);
% exportgraphics(gcf,strcat(path,image_name,'/prints/','cuprite_compare_all_noise0_R30_4algs','.pdf'))

% Update metrics
H_FCLSU = bundle{1,k1}*A_FCLSU_all{1,k1};
H_collaborative = bundle{1,k2}*A_collaborative_all{1,k2};
H_group = bundle{1,k3}*A_group_all{1,k3};
H_elitist = bundle{1,k4}*A_elitist_all{1,k4}; 
H_fractional = bundle{1,k5}*A_fractional_all{1,k5};
H_SUnCNN = bundle{1,k7}*A_SUnCNN_all{1,k7};
H_MUA = bundle{r,k6}*A_MUA_all{r,k6}; 

SRE_FCLSU = 20*log10(norm(H_true,'fro')/norm(H_FCLSU - H_true,'fro'));
SRE_collaborative =  20*log10(norm(H_true,'fro')/norm(H_collaborative - H_true,'fro'));
SRE_group = 20*log10(norm(H_true,'fro')/norm(H_group - H_true,'fro'));
SRE_elitist = 20*log10(norm(H_true,'fro')/norm(H_elitist - H_true,'fro'));
SRE_fractional = 20*log10(norm(H_true,'fro')/norm(H_fractional - H_true,'fro'));
SRE_SUnCNN = 20*log10(norm(H_true,'fro')/norm(H_SUnCNN - H_true,'fro'));
SRE_MUA = 20*log10(norm(H_true,'fro')/norm(H_MUA - H_true,'fro'));

% RMSE
RMSE_FCLSU = sqrt(1/L*sum((H_FCLSU-X).^2,1));
RMSE_collaborative =  sqrt(1/L*sum((H_collaborative-X).^2,1));
RMSE_group = sqrt(1/L*sum((H_group-X).^2,1));
RMSE_elitist = sqrt(1/L*sum((H_elitist-X).^2,1));
RMSE_fractional = sqrt(1/L*sum((H_fractional-X).^2,1));
RMSE_SUnCNN = sqrt(1/L*sum((H_SUnCNN-X).^2,1));
RMSE_MUA = sqrt(1/L*sum((H_MUA-X).^2,1));

% SAM
SAM_true = zeros(N,1);
SAM_FCLSU = zeros(N,1);
SAM_collaborative = zeros(N,1);
SAM_group = zeros(N,1);
SAM_elitist = zeros(N,1);
SAM_fractional = zeros(N,1);
SAM_SUnCNN = zeros(N,1);
SAM_MUA = zeros(N,1);

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

for k = 1:N
    SAM_SUnCNN(k) = 180/pi*real(acos((X(:,k)'*H_SUnCNN(:,k))...
        /(norm(X(:,k))*norm(H_SUnCNN(:,k)))));
end

for k = 1:N
    SAM_MUA(k) = 180/pi*real(acos((X(:,k)'*H_MUA(:,k))...
        /(norm(X(:,k))*norm(H_MUA(:,k)))));
end

%% Display RMSEs, SAMs, SREs

fprintf('\n--- SAMs \n\n')
fprintf('FCLSU = %g \n', mean(SAM_FCLSU(:)))
fprintf('collaborative = %g \n', mean(SAM_collaborative(:)))
fprintf('group = %g \n', mean(SAM_group(:)))
fprintf('elitist = %g \n', mean(SAM_elitist(:)))
fprintf('fractional = %g \n', mean(SAM_fractional(:)))
fprintf('SUnCNN = %g \n', mean(SAM_SUnCNN(:)))
fprintf('GMBUA = %g \n', SAM_MUA(MUA_median_idx))

fprintf('\n--- RMSEs \n\n')
fprintf('FCLSU = %g \n', mean(RMSE_FCLSU(:)))
fprintf('collaborative = %g \n', mean(RMSE_collaborative(:)))
fprintf('group = %g \n', mean(RMSE_group(:)))
fprintf('elitist = %g \n', mean(RMSE_elitist(:)))
fprintf('fractional = %g \n', mean(RMSE_fractional(:)))
fprintf('SUnCNN = %g \n', mean(RMSE_SUnCNN(:)))
fprintf('GMBUA = %g \n', RMSE_MUA(MUA_median_idx))

fprintf('\n--- SREs - Observed Image \n\n')
fprintf('FCLSU = %g \n', SRE_FCLSU_median)
fprintf('collaborative = %g \n', SRE_collaborative_median)
fprintf('group = %g \n', SRE_group_median)
fprintf('elitist = %g \n', SRE_elitist_median)
fprintf('fractional = %g \n', SRE_fractional_median)
fprintf('SUnCNN = %g \n', SRE_SUnCNN_median)
fprintf('GMBUA = %g \n', SRE_MUA_median)

fprintf('\n--- Exec. - Time \n\n')
fprintf('SUnCNN = %g \n', median(time_SUnCNN))

% %% Global abundance maps -- All
% 
% for r = 1:30
%     A_FCLSU_im = reshape(A_FCLSU_global_all{1,r}',m,n,P);
%     A_collaborative_im =reshape(A_collaborative_global_all{1,r}',m,n,P);
%     A_group_im = reshape(A_group_global_all{1,r}',m,n,P);
%     A_elitist_im = reshape(A_elitist_global_all{1,r}',m,n,P);
%     A_fractional_im = reshape(A_fractional_global_all{1,r}',m,n,P);
%     A_MUA_im = reshape(A_MUA_global_all{r,r}',m,n,P);
% 
%     vertical_spacing = 0.002;
%     horizontal_spacing = 0.002;
%     font_s = 13;
%     figure, set(gcf,'color', 'white')
%     for p = 1:P
%         subaxis(6,P,p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
%         imshow(A_FCLSU_im(:,:,p),[],'colormap', jet)
%         set(gca,'clim',[0,1]), axis square
%         if p == 1
%             ylabel('FCLSU','fontname','times','fontsize',font_s)
%         end
%         colormap jet
% 
%         subaxis(6,P,P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
%         imshow(A_collaborative_im(:,:,p),[],'colormap', jet)
%         set(gca,'clim',[0,1]), axis square
%         if p == 1
%             ylabel('Collab.','fontname','times','fontsize',font_s)
%         end
%         colormap jet
% 
%         subaxis(6,P,2*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
%         imshow(A_group_im(:,:,p),[],'colormap', jet)
%         set(gca,'clim',[0,1]), axis square
%         if p == 1
%             ylabel('Group','fontname','times','fontsize',font_s)
%         end
%         colormap jet
% 
%         subaxis(6,P,3*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
%         imshow(A_elitist_im(:,:,p),[],'colormap', jet)
%         set(gca,'clim',[0,1]), axis square
%         if p == 1
%             ylabel('Elitist','fontname','times','fontsize',font_s)
%         end
%         colormap jet
% 
%         subaxis(6,P,4*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
%         imshow(A_fractional_im(:,:,p),[],'colormap', jet)
%         set(gca,'clim',[0,1]), axis square
%         if p == 1
%             ylabel('Fractional','fontname','times','fontsize',font_s)
%         end
%         colormap jet
% 
%         subaxis(6,P,5*P+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
%         imshow(A_MUA_im(:,:,p),[],'colormap', jet)
%         set(gca,'clim',[0,1]), axis square
%         if p == 1
%             ylabel('GMBUA','fontname','times','fontsize',font_s)
%             xlabel('Em 1','fontname','times','fontsize',font_s)
%         elseif p == 2
%             xlabel('Em 2','fontname','times','fontsize',font_s)
%         elseif p == 3
%             xlabel('Em 3','fontname','times','fontsize',font_s)
%         elseif p == 4
%             xlabel('Em 4','fontname','times','fontsize',font_s)
%         elseif p == 5
%             xlabel('Em 5','fontname','times','fontsize',font_s)
%         elseif p == 6
%             xlabel('Em 6','fontname','times','fontsize',font_s)
%         elseif p == 7
%             xlabel('Em 7','fontname','times','fontsize',font_s)
%         elseif p == 8
%             xlabel('Em 8','fontname','times','fontsize',font_s)
%         elseif p == 9
%             xlabel('Em 9','fontname','times','fontsize',font_s)
%         elseif p == 10
%             xlabel('Em 10','fontname','times','fontsize',font_s)
%         elseif p == 11
%             xlabel('Em 11','fontname','times','fontsize',font_s)
%         elseif p == 12
%             xlabel('Em 12','fontname','times','fontsize',font_s)
%         elseif p == 13
%             xlabel('Em 13','fontname','times','fontsize',font_s)
%         elseif p == 14
%             xlabel('Em 14','fontname','times','fontsize',font_s)
%         end
%     end
% end
% 
% %% Global abundance maps -- Manual select -- All
% 
% for r = 1:30
% 
%     A_FCLSU_im = reshape(A_FCLSU_global_all{1,r}',m,n,P);
%     A_collaborative_im =reshape(A_collaborative_global_all{1,r}',m,n,P);
%     A_group_im = reshape(A_group_global_all{1,r}',m,n,P);
%     A_elitist_im = reshape(A_elitist_global_all{1,r}',m,n,P);
%     A_fractional_im = reshape(A_fractional_global_all{1,r}',m,n,P);
%     A_MUA_im = reshape(A_MUA_global_all{r,r}',m,n,P);
% 
%     vertical_spacing = 0.002;
%     horizontal_spacing = 0.002;
%     font_s = 13;
%     figure, set(gcf,'color', 'white')
%     for p = 1:length(ems)
%         subaxis(6,length(ems),p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
%         imshow(A_FCLSU_im(:,:,ems(p)),[],'colormap', jet)
%         set(gca,'clim',[0,1]), axis square
%         if p == 1
%             ylabel('FCLSU','fontname','times','fontsize',font_s)
%         end
%         colormap jet
% 
%         subaxis(6,length(ems),length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
%         imshow(A_collaborative_im(:,:,ems(p)),[],'colormap', jet)
%         set(gca,'clim',[0,1]), axis square
%         if p == 1
%             ylabel('Collab.','fontname','times','fontsize',font_s)
%         end
%         colormap jet
% 
%         subaxis(6,length(ems),2*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
%         imshow(A_group_im(:,:,ems(p)),[],'colormap', jet)
%         set(gca,'clim',[0,1]), axis square
%         if p == 1
%             ylabel('Group','fontname','times','fontsize',font_s)
%         end
%         colormap jet
% 
%         subaxis(6,length(ems),3*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
%         imshow(A_elitist_im(:,:,ems(p)),[],'colormap', jet)
%         set(gca,'clim',[0,1]), axis square
%         if p == 1
%             ylabel('Elitist','fontname','times','fontsize',font_s)
%         end
%         colormap jet
% 
%         subaxis(6,length(ems),4*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
%         imshow(A_fractional_im(:,:,ems(p)),[],'colormap', jet)
%         set(gca,'clim',[0,1]), axis square
%         if p == 1
%             ylabel('Fractional','fontname','times','fontsize',font_s)
%         end
%         colormap jet
% 
%         subaxis(6,length(ems),5*length(ems)+p,'SpacingVertical',vertical_spacing,'SpacingHorizontal',horizontal_spacing)
%         imshow(A_MUA_im(:,:,ems(p)),[],'colormap', jet)
%         set(gca,'clim',[0,1]), axis square, axis square
%         if p == 1
%         ylabel('GMBUA','fontname','times','fontsize',font_s)
%         xlabel('Alunite','fontname','times','fontsize',font_s)
%         elseif p == 2
%             xlabel('Muscovite','fontname','times','fontsize',font_s)
%         elseif p == 3
%             xlabel('Kaolinite','fontname','times','fontsize',font_s)
%         elseif p == 4
%             xlabel('Sphene','fontname','times','fontsize',font_s)
%         elseif p == 5
%             xlabel('Buddingtonite','fontname','times','fontsize',font_s)
%         elseif p == 6
%             xlabel('Chalcedony','fontname','times','fontsize',font_s)
%         end
%     end
% end

