%% Notes
% Hue = (650 - wavelength)*240/(650-475);
% L = 620 - 170 / 270 * H

%% Get image files
clear all;
close all;

imgDir = uigetdir;
imgDir = [imgDir '/'];
imgList=[dir([imgDir '*.jpg']);dir([imgDir '*.png']);dir([imgDir '*.tif'])];
n_img=length(imgList);
% n_group = 2;
% n_indiv = 5;
 n_group = 4;
 n_indiv = 9;
n_color = 256;
cmap = hsv(n_color);

%%

for i = 1:n_img
%     figure;
    rgbImage = imread([imgDir imgList(i).name]);
    [rows, columns, numberOfColorBands] = size(rgbImage);

    mask = rgbImage(:,:,1) >25 & rgbImage(:,:,2) < 240 & rgbImage(:,:,3) >25;
    maskedRgbImage = bsxfun(@times, rgbImage, cast(mask, class(rgbImage)));
    
    % Get hue channel
    hsvImage = rgb2hsv(maskedRgbImage);
    h = hsvImage(:,:,1);
    h = h(h~=0);
    [hueHist(i,:), grayLevels(i,:)] = hist(h(:), n_color);
    hueHistNorm(i,:) = hueHist(i,:)/norm(hueHist(i,:));
end

for i = 1:n_group
groupHist(i,:) = sum(hueHist((i-1)*n_indiv+1:(i-1)*n_indiv+n_indiv,:));
groupHistNorm(i,:) = groupHist(i,:)/norm(groupHist(i,:));
end

% for i = 1:n_group
%     figure;
%     hold on;
%     for j = 1:n_color
%     patch([j-.5,j+.5,j+.5,j-.5],[0,0,groupHist(i,j),groupHist(i,j)],cmap(j,:),'FaceAlpha', 1);
%     end
%     hold off;
%     grid on;
%     title(['Histogram of Hue for T' num2str(i-1)], 'FontSize', 20);
%     set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% end


%% Overview of Group Color Histogram
i_hMean = 1;
hMean = [];
hErr = [];
for i = 1:n_indiv:n_group*n_indiv
    hMean(i_hMean,:) = mean(hueHistNorm(i:i+n_indiv-1,:));
    hErr(i_hMean,:) = std(hueHistNorm(i:i+n_indiv-1,:))/sqrt(n_indiv);
    i_hMean = i_hMean+1;
end

i_c = find(sum(hMean,1) >0.0005);
i_c1 = find(hMean(1,:) >0.0005);
i_c2 = find(hMean(2,:) >0.0005);
i_c3 = find(hMean(3,:) >0.0005);
i_c4 = find(hMean(4,:) >0.0005);

hold on;
for i = 1:length(i_c)
    patch([i-.5,i+.5,i+.5,i-.5],[0,0,-.01,-.01],cmap(i_c(i),:),'FaceAlpha', 1);
end
p1 = plot(hMean(1,i_c),'k^-','MarkerSize',8,'MarkerFaceColor','r');
p2 = plot(hMean(2,i_c),'ko-','MarkerSize',8,'MarkerFaceColor','g');
p3 = plot(hMean(3,i_c),'ks-','MarkerSize',8,'MarkerFaceColor','b');
p4 = plot(hMean(4,i_c),'kp-','MarkerSize',8,'MarkerFaceColor','y');
hold off;

% lgd = legend([p1,p2],...
%     'Yellow Group','Blue Group',...
%     'Location','Best');
lgd = legend([p1,p2],...
    'IRL','PBS',...
    'Location','Best');
xlabel('Colors'), ylabel('Count');
ylim([-0.01,Inf]);
title('Histogram of Before and After Treatment for IRL and PBS');
% close;
%% Which color passes ANOVA?
hueHistNorm3D = reshape(hueHistNorm,n_indiv,n_group,n_color);

i_cSig = [];
for j = 1:n_color
    p = anova1(hueHistNorm3D(:,:,j),[],'off');
    %p = kruskalwallis(hueHistNorm3D(:,:,j),[],'off');
    if p < 0.00001 && max(mean(hueHistNorm3D(:,:,j))) > 0.05
        i_cSig(end+1) = j;
    end
end

i_cSig12 = [];
for j = i_cSig
   [h,p12] =ttest2(hueHistNorm3D(:,2,j), hueHistNorm3D(:,3,j));
   if p12 < 0.001
       i_cSig12(end+1) = j;
   end
end

figure;
hold on;
for i = 1:length(i_cSig12)
    patch([i-.5,i+.5,i+.5,i-.5],[0,0,-.01,-.01],cmap(i_cSig12(i),:),'FaceAlpha', 1);
end
p2Sig = errorbar(hMean(2,i_cSig12),hErr(2,i_cSig12));
p3Sig = errorbar(hMean(3,i_cSig12),hErr(3,i_cSig12));
p2Sig = area(hMean(2,i_cSig12),'FaceColor','b','FaceAlpha',0.4);
p3Sig = area(hMean(3,i_cSig12),'FaceColor','y','FaceAlpha',0.4);
hold off;
lgd = legend([p2Sig,p3Sig],...
    'T1','T2',...
    'Location','Best');
set(lgd,'FontSize',30);
xlabel('Colors', 'FontSize', 16), ylabel('Proportion', 'FontSize', 16);
ylim([-0.01,Inf]);
title('Pair-wise T-Test between T1 and T2 (p < 0.01)','fontsize',20);

i_cSig23 = [];
for j = i_cSig
   [h,p23] =ttest2(hueHistNorm3D(:,3,j), hueHistNorm3D(:,4,j));
   if p23 < 0.001
       i_cSig23(end+1) = j;
   end
end

figure;
hold on;
for i = 1:length(i_cSig23)
    patch([i-.5,i+.5,i+.5,i-.5],[0,0,-.01,-.01],cmap(i_cSig23(i),:),'FaceAlpha', 1);
end
p3Sig = errorbar(hMean(3,i_cSig23),hErr(3,i_cSig23));
p4Sig = errorbar(hMean(4,i_cSig23),hErr(4,i_cSig23));
p3Sig = area(hMean(3,i_cSig23),'FaceColor','y','FaceAlpha',0.4);
p4Sig = area(hMean(4,i_cSig23),'FaceColor','b','FaceAlpha',0.4);
hold off;
lgd = legend([p3Sig,p4Sig],...
    'T2','T3',...
    'Location','Best');
set(lgd,'FontSize',30);
xlabel('Colors', 'FontSize', 16), ylabel('Proportion', 'FontSize', 16);
ylim([-0.01,Inf]);
title('Pair-wise T-Test between T2 and T3 (p < 0.01)','fontsize',20);

i_cSig13 = [];
for j = i_cSig
   [h,p13] =ttest2(hueHistNorm3D(:,2,j), hueHistNorm3D(:,4,j));
   if p13 < 0.001
       i_cSig13(end+1) = j;
   end
end

figure;
hold on;
for i = 1:length(i_cSig13)
    patch([i-.5,i+.5,i+.5,i-.5],[0,0,-.01,-.01],cmap(i_cSig13(i),:),'FaceAlpha', 1);
end
p2Sig = errorbar(hMean(2,i_cSig13),hErr(2,i_cSig13));
p4Sig = errorbar(hMean(4,i_cSig13),hErr(4,i_cSig13));
p2Sig = area(hMean(2,i_cSig13),'FaceColor','b','FaceAlpha',0.4);
p4Sig = area(hMean(4,i_cSig13),'FaceColor','r','FaceAlpha',0.4);
hold off;
lgd = legend([p2Sig,p4Sig],...
    'T1','T3',...
    'Location','Best');
set(lgd,'FontSize',30);
xlabel('Colors', 'FontSize', 16), ylabel('Proportion', 'FontSize', 16);
ylim([-0.01,Inf]);
title('Pair-wise T-Test between T1 and T3 (p < 0.01)','fontsize',20);


figure;
hold on;
for i = 1:length(i_cSig)
    patch([i-.5,i+.5,i+.5,i-.5],[0,0,-.01,-.01],cmap(i_cSig(i),:),'FaceAlpha', 1);
end

p2Sig = errorbar(hMean(2,i_cSig),hErr(2,i_cSig));
p3Sig = errorbar(hMean(3,i_cSig),hErr(3,i_cSig));
p4Sig = errorbar(hMean(4,i_cSig),hErr(4,i_cSig));
p2Sig = area(hMean(2,i_cSig),'FaceColor','b','FaceAlpha',0.4);
p3Sig = area(hMean(3,i_cSig),'FaceColor','y','FaceAlpha',0.4);
p4Sig = area(hMean(4,i_cSig),'FaceColor','r','FaceAlpha',0.4);

hold off;
lgd = legend([p2Sig,p3Sig,p4Sig],...
    'T1','T2','T3',...
    'Location','Best');
set(lgd,'FontSize',30);
xlabel('Colors', 'FontSize', 16), ylabel('Proportion', 'FontSize', 16);
ylim([-0.01,Inf]);
title('ANOVA among T1, T2 and T3 (p < 0.00001)','fontsize',20);


%% Chi-Square Histogram Distance
% dist_func=@chi_square_statistics_fast;
% D=pdist2(hueHistNorm,hueHistNorm,dist_func);
% hm = HeatMap(D);
% hm.addTitle('Chi-Square Distance Heat Map');
% hm.plot;
figure;
% Eigenspace Analysis
cSig = cmap(i_cSig,:);
c1 = cmap(i_c1,:);
c2 = cmap(i_c2,:);
c3 = cmap(i_c3,:);
c4 = cmap(i_c4,:);
cSigM = mean(cSig);
cSigV = cov(cSig);
[V,D] = eig(cSigV);
S = zeros(3,3);
S(:,1) = V(:,3)*D(3,3);
S(:,2) = V(:,2)*D(2,2);
S(:,3) = V(:,1)*D(1,1);

S = S/norm(S);

scatter3(cSig(:,1),cSig(:,2),cSig(:,3),200,cSig,'filled');
xlabel('red'), ylabel('green'), zlabel('blue');
title('Most Significantly Different Colors with Eigenvectors');
axis equal;
view(19,31);
hold on;
quiver3(repmat(cSigM(1),1,3),repmat(cSigM(2),1,3),repmat(cSigM(3),1,3),...
    S(1,:),S(2,:),S(3,:));
hold off;
grid on;

% plot mean color point clouds for T0, T1, T2, T3
figure;
hold on;
scatter3(c1(:,1) + (rand(length(c1(:,1)),1)-0.5)*0.05,...
    c1(:,2) + (rand(length(c1(:,1)),1)-0.5)*0.05,...
    c1(:,3) + + (rand(length(c1(:,1)),1)-0.5)*0.05,...
    100,'r','filled','^');
scatter3(c2(:,1) + (rand(length(c2(:,1)),1)-0.5)*0.05,...
    c2(:,2) + (rand(length(c2(:,1)),1)-0.5)*0.05,...
    c2(:,3) + (rand(length(c2(:,1)),1)-0.5)*0.05,...
    100,'g','filled','s');
scatter3(c3(:,1) + (rand(length(c3(:,1)),1)-0.5)*0.05,...
    c3(:,2) + (rand(length(c3(:,1)),1)-0.5)*0.05,...
    c3(:,3) + (rand(length(c3(:,1)),1)-0.5)*0.05,...
    100,'b','filled','^');
scatter3(c4(:,1) + (rand(length(c4(:,1)),1)-0.5)*0.05,...
    c4(:,2) + (rand(length(c4(:,1)),1)-0.5)*0.05,...
    c4(:,3) + (rand(length(c4(:,1)),1)-0.5)*0.05,...
    100,'y','filled','^');
legend('T0','T1','T2','T3','Location','Best');
% legend('Yellow Group','Blue Grouop','Location','Best');
quiver3(repmat(cSigM(1),1,3),repmat(cSigM(2),1,3),repmat(cSigM(3),1,3),...
    S(1,:),S(2,:),S(3,:));
hold off;
xlabel('red'), ylabel('green'), zlabel('blue');
title('Mean Per Group Colors with Eigenvectors');
axis equal;
view(17,18);
grid on;
figure;

% Principal Component Color Complements
for i = 1:3
    ind = dsearchn(cSig,S(:,i)'*10);
    eigcolors1(i,:)=cSig(ind,:);
    ind = dsearchn(cSig,-S(:,i)'*10);
    eigcolors2(i,:)=cSig(ind,:);
end

for i = 1:3
%     patch([i-1 i i i-1], [0 0 0.5 0.5], (repmat(0.5,3,1)+S(:,i)/norm(S(:,i))/2)')
%     patch([i-1 i i i-1], [0.5 0.5 1 1], (repmat(0.5,3,1)-S(:,i)/norm(S(:,i))/2)')
    patch([i-1 i i i-1], [0 0 0.5 0.5], eigcolors1(i,:))
    patch([i-1 i i i-1], [0.5 0.5 1 1], eigcolors2(i,:))
end

