% select folder that contains images to analyze
PathName = uigetdir;
% might have to switch from / to \ for Windows
PathName = [PathName '/'];
FileList=[dir([PathName '*.jpg']);dir([PathName '*.png']);dir([PathName '*.tiff'])];

%% histogram match within groups
% for i = 1:length(FileList)
%     if mod(i,9) == 1
%     ref = imread([PathName FileList(i).name]);
%     imwrite(ref, [PathName 'new_' FileList(i).name]);
%     else
%     tempim = imread([PathName FileList(i).name]);
%     im = imhistmatch(tempim,ref);
%     imwrite(im, [PathName 'new_' FileList(i).name]);
%     end
% end

%% luminosity match across groups
meanlumin = 0;
for i = 1:length(FileList)
    rgb = imread([PathName FileList(i).name]);
    lab = rgb2lab(rgb);
    meanlumin = meanlumin + mean2(lab(:,:,1));
    disp(['Progress: ' num2str(i) '/' num2str(length(FileList))]);
end
meanlumin = meanlumin / length(FileList);

for i = 1:length(FileList)
    rgb = imread([PathName FileList(i).name]);
    lab = rgb2lab(rgb);
    lab(:,:,1) = lab(:,:,1) + ( meanlumin - mean2(lab(:,:,1)) );
    rgb = lab2rgb(lab);
    imwrite(rgb, [PathName 'final_' FileList(i).name]);
    disp(['Progress: ' num2str(i) '/' num2str(length(FileList))]);
end

%%

% number of maximum color differences to show
num_of_max = 10;
bin_size = 51;

col_val_arr = [];
col_freq_arr = [];
col_all_possible = [];
col_all_R = [];
col_all_G = [];
col_all_B = [];

for i = 0:bin_size:255
    for j = 0:bin_size:255
        for k = 0:bin_size:255
            col_all_possible(end+1,1) = bin_size*floor(k/bin_size)+256*(bin_size*floor(j/bin_size))+256*256*(bin_size*floor(i/bin_size));
            col_all_R(end+1,1) = k;
            col_all_G(end+1,1) = j;
            col_all_B(end+1,1) = i;
        end
    end
end

% compensate for the number of bins used in histogram
col_all_possible = col_all_possible + 1;
col_all_possible(1) = col_all_possible(1) - 1;
T_color = table(col_all_possible(2:end),col_all_R(2:end),...
   col_all_G(2:end),col_all_B(2:end));
T_color.Properties.VariableNames = {'col_category','R','G','B'};

for i = 1:length(FileList)
im_rgb =double(imread([PathName FileList(i).name]));
[~,name,~] = fileparts([PathName FileList(i).name]);
col_count=bin_size*floor(im_rgb(:,:,1)/bin_size)+256*(bin_size*floor(im_rgb(:,:,2)/bin_size))+256*256*(bin_size*floor(im_rgb(:,:,3)/bin_size));
N = histcounts(col_count,col_all_possible);
T = table(N');
T.Properties.VariableNames = {name};
T_color = [T_color, T];
disp(['Progress: ' num2str(i)]);
end
[~, index] = max(T_color.T0_fish1_body);
T_color(index,:) = [];

writetable(T_color,[PathName 'color_analysis.csv']);

num_of_max = 20;
v1 = T_color{:,5+9};
v2 = T_color{:,5+18};
[sortedX,sortingIndices] = sort(abs(v1-v2)./abs(v1+v2),'descend');
sortingIndices = sortingIndices(~isnan(sortedX) & sortedX ~= 1);
col_diff = T_color.col_category(sortingIndices(1:num_of_max));
col_diff = col_diff-1;
col_val_rgb=([rem(col_diff,256) floor(rem(col_diff,256*256)/256) floor(col_diff/(256*256))]/255);

for i = 1:num_of_max
    patch([i-1,i,i,i-1],[0,0,1,1],col_val_rgb(i,:));
end

%% Check color by value
whatcolor = 65280;
col_val_rgb=([rem(whatcolor,256) floor(rem(whatcolor,256*256)/256) floor(whatcolor/(256*256))]/255);
patch([0,1,1,0],[0,0,1,1],col_val_rgb);

%%
% Between T0,1,2,3 group color mean profile analysis
% for each bin color, make table of columns with T0,1,2,3 and rows with
% each fish's color count profile in that color bin.
% perform kruskal-wallis test to compare each color profile
% this should result in 63 color profiles
% find significant profiles and get their colors