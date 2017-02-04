% [1] B. Schauerte, G. A. Fink, "Web-based Learning of Naturalized Color 
%     Models for Human-Machine Interaction". In Proceedings of the 12th 
%     International Conference on Digital Image Computing: Techniques and 
%     Applications (DICTA), IEEE, Sydney, Australia, December 1-3, 2010. 
% [2] B. Schauerte, R. Stiefelhagen, "Learning Robust Color Name Models 
%     from Web Images". In Proceedings of the 21st International Conference
%     on Pattern Recognition (ICPR), Tsukuba, Japan, November 11-15, 2012

PathName = uigetdir;
PathName = [PathName '/'];
FileList=[dir([PathName '*.jpg']);dir([PathName '*.png']);dir([PathName '*.tiff'])];
n_samples=length(FileList);

%% Calculate color image histograms
n_bins=10;
edges=(0:n_bins)/n_bins;
edges(end) = edges(end) + 0.0000001;
histograms=zeros(n_samples,n_bins*n_bins*n_bins);
histograms_raw=zeros(n_samples,n_bins*n_bins*n_bins);

for i=1:n_samples
  I=imread([PathName FileList(i).name]);
  I=im2double(I); 
  [~,r_bins] = histc(reshape(I(:,:,1),1,[]),edges);
  [~,g_bins] = histc(reshape(I(:,:,2),1,[]),edges); 
  [~,b_bins] = histc(reshape(I(:,:,3),1,[]),edges); 
  histogram=zeros(n_bins,n_bins,n_bins); 
  
  for j=1:numel(r_bins)
    histogram(r_bins(j),g_bins(j),b_bins(j)) = histogram(r_bins(j),g_bins(j),b_bins(j)) + 1;
  end 
  histograms(i,:) = reshape(histogram,1,[]) / sum(histogram(:));
  histograms_raw(i,:) = reshape(histogram,1,[]);
  
  disp(['Progress: ' num2str(i) '/' num2str(n_samples)]);
end

% compensate for fish size by multiplying the ratio of background colors
[~,ind] = max(histograms(1,:));
scale = histograms(:,ind) / max(histograms(:,ind));
histograms = histograms .* scale;
histograms(:,ind) = [];
%% Overview of color histogram
ind = 1;
totalhist = [];
for i = [1, 10, 19, 28]
    totalhist(ind,:) = mean(histograms(i:i+9-1,:));
    ind = ind+1;
end

meanedges = edges(2:end) - diff(edges)/2;
ind = 1;
for b = 1:n_bins
    for g = 1:n_bins
        for r = 1:n_bins
            col_val_rgb(ind,:) = [meanedges(r), meanedges(g), meanedges(b)];
            ind = ind + 1;
        end
    end
end

hold on;
for i = 1:length(col_val_rgb)
    patch([i-.5,i+.5,i+.5,i-.5],[0,0,.02,.02],col_val_rgb(i,:),'FaceAlpha', 1);
end

plot(totalhist(1,:),'w^-');
plot(totalhist(2,:),'wo-');
plot(totalhist(3,:),'ws-');
plot(totalhist(4,:),'w*-');
hold off;

%% Which color passes ANOVA?
for j = 1:size(histograms,2)
    ind = 1;
for i = [1, 10, 19, 28]
    bygroup(:,ind,j) = histograms(i:i+9-1,j);
    ind = ind+1;
end
end

ind = [];
for j = 1:size(histograms,2)
    p = anova1(bygroup(:,:,j),[],'off');
    if p < 0.001 && max(mean(bygroup(:,:,j))) > 0.0005
        ind(end+1) = j;
    end
end

for i = 1:length(ind)
    patch([i-.5,i+.5,i+.5,i-.5],[0,0,.02,.02],col_val_rgb(ind(i),:),'FaceAlpha', 1);
end


%%
dist_func=@chi_square_statistics_fast;
%dist_func=@chi_square_statistics;
D=pdist2(histograms,histograms,dist_func);
HeatMap(D);