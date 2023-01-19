close all;
clear all;

% Add helper functions to path
addpath('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/connectomes-data/Toolboxes/2016_01_16_BCT/');
addpath('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/connectomes-data/Toolboxes/NIfTI_20140122/');
addpath('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/connectomes-data/Toolboxes/covshrink-kpm/');
addpath('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/connectomes-data/Toolboxes/hex_and_rgb_v1/');

FA_threshold = [0.1:0.1:0.8];
% Load nii files
filename = '/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/connectomes-data/Task1Data/tractor/diffusion/dti_FA.nii.gz';
result = load_nii(filename);

% Get the fa map from the results object
fa_images = result.img;
fa_images(isnan(fa_images)) = 0;

% Plot original image
slice_num = 15;
subplot(3,3,1);
imshow(fa_images(:,:,slice_num))
title('Original')
%%
parcellation_func = load_nii('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/connectomes-data/Task1Data/tractor/diffusion/parcellation.nii.gz');
result = load_nii('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/connectomes-data/Task1Data/tractor/diffusion/dti_FA.nii.gz');
color_region_file = fopen('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/connectomes-data/Task1Data/tractor/diffusion/parcellation.lut');
region = textscan(color_region_file, '%d %s %s %s %s %s %s', 'Headerlines', 3); % start from foutrh line

fa_images = result.img;
fa_images(isnan(fa_images)) = 0;

%% doing this for the original FA threshood = 0 
parce_images = parcellation_func.img;
overlay_img = zeros(96, 96, 3);
slice_num = 40;

% for i=1:96
%     for j=1:96
%         if parce_images(i,j,slice_num) ~= 0
%             int_to_index = find([region{1,1}] == parce_images(i,j,slice_num));
%             index_to_hex = region{1,4}{int_to_index};
%             hex_to_rgb = hex2rgb(index_to_hex(2:end-1));
%             overlay_img(i,j, :) = hex_to_rgb;
%         else
%             overlay_img(i,j, :) = fa_images(i,j,slice_num);
%         end 
%     end 
% end 
% imshow(overlay_img)

% slice_num = 15;
% subplot(3,3,1);
% temp = parce_images(:,:,slice_num);
% imshow(parce_images(:,:,slice_num))
% title('Original')
% 
% % Plot original image
% slice_num = 15;
% subplot(3,3,1);
imshow(fa_images(:,:,slice_num))
% title('Original')
  
%% question 2
count_subvoxel = zeros([1,8]);
for i = 1:8
    % Threshold fa map
    threshold = 0.1 * i;
    fa_images_threshold = fa_images;
    fa_images_threshold(fa_images_threshold < threshold) = 0;
    
    % Plot
    subplot(3,3,1+i);
    imshow(fa_images_threshold(:,:,slice_num));

    slice_img = fa_images_threshold(:,:,slice_num);
    count_subvoxel(i) = length(slice_img(slice_img>0));
    title(threshold);
end

figure
hold on 
bar(FA_threshold,count_subvoxel)
xlabel('FA threshold')
ylabel('The number of voxels thresholded out')

%% question 3
connectome_file = {'FA.1_graph.csv'; 'FA.2_graph.csv'; 'FA.3_graph.csv';...
            'FA.4_graph.csv'; 'FA.5_graph.csv'; 'FA.6_graph.csv'; ...
            'FA.7_graph.csv'; 'FA.8_graph.csv'};

density = zeros([1,8]);
char_path = zeros([1,8]);
efficiency = zeros([1,8]);
mean_cluster_coeff = zeros([1,8]);
for i = 1:8
    fa_file_path = ['/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/connectomes-data/Task1Data/connectomes/', connectome_file{i}];
    % skip the first three lines 
    struct_connect{i} = csvread(fa_file_path,3,0);
    [density(i),char_path(i),efficiency(i),mean_cluster_coeff(i)] = graph_metrics(struct_connect{i});
end
display_graphs(density,char_path,efficiency,mean_cluster_coeff, ...
    'FA Threshold','When FA threshold is changing:');

%% question 4
parcellation_func = load_nii('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/connectomes-data/Task1Data/tractor/functional/parcellation.nii.gz');
data_func = load_nii('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/connectomes-data/Task1Data/tractor/functional/data.nii.gz');
% imshow(parcellation_func(:,:,1))

% Load cortical region information
file_id = fopen('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/connectomes-data/Task1Data/tractor/functional/parcellation.lut');
region = textscan(file_id, '%d %s %s %s %s %s %s', 'Headerlines', 3); % start from foutrh line
region_ids = region{1}; 
[num_regions,~] = size(region_ids);

% now we are going to compute the average time series for each cortical region
% the size of rs-fMRI is 64×64×30×15, and the fourth dimension is time 
avg_time_series = zeros(num_regions, 15);
for i = 1:num_regions
    id = region_ids(i);

    % use mask to select voxels of this cortical region
    cortical_region = (parcellation_func.img == id); 
    voxels = data_func.img(cortical_region(:,:,:,ones(1,15)));
    [size_voxels,~] = size(voxels); 
    region_voxels = reshape(voxels, [size_voxels/15, 15]); 
    
    % average time series of this cortical region's all voxel
    avg_time_series(i,:) = mean(region_voxels, 1);
end

% corshrink
for i = 1:8
    corshrink_matrix{i} = corshrink(avg_time_series',FA_threshold(i));
end

% resulting matrics at a correlation of 0.1, binarise it
for t = 1:8
    corshrink_new = corshrink_matrix{t};
    % binarise it
    corshrink_new(corshrink_new < 0.1) = 0;
    corshrink_new(corshrink_new >= 0.1) = 1;
    binary_corshrink{t} = corshrink_new;
end

%  graph metrics 
cor_density = zeros([1,8]);
cor_char_path = zeros([1,8]);
cor_efficiency = zeros([1,8]);
cor_mean_cluster_coeff = zeros([1,8]);
for i = 1:8
    [cor_density(i),cor_char_path(i),cor_efficiency(i),cor_mean_cluster_coeff(i)] = ...
        graph_metrics(binary_corshrink{i});
end
display_graphs(cor_density,cor_char_path,cor_efficiency,cor_mean_cluster_coeff, ...
    'lambda','When lambda is changing:');

% discarding negative correlations with absolute value greater than 0.1
for t = 1:8
    corshrink_new = corshrink_matrix{t};
    corshrink_new = (corshrink_new > -0.1);
    discard_neg_corshrink{t} = corshrink_new;
end

discard_neg_density = zeros([1,8]);
discard_neg_char_path = zeros([1,8]);
discard_neg_efficiency = zeros([1,8]);
discard_neg_mean_cluster_coeff = zeros([1,8]);
for i = 1:8
    [discard_neg_density(i),discard_neg_char_path(i),discard_neg_efficiency(i),...
        discard_neg_mean_cluster_coeff(i)] = graph_metrics(discard_neg_corshrink{i});
end

show_difference_graphs(cor_density,discard_neg_density,'edge density');
show_difference_graphs(cor_char_path,discard_neg_char_path,'mean shortest path');
show_difference_graphs(cor_efficiency,discard_neg_efficiency,'efficiency');
show_difference_graphs(cor_mean_cluster_coeff,discard_neg_mean_cluster_coeff,'mean cluster coeff');

% function definition
% mannually setting lambda
function [Rhat] = corshrink(x,input_lambda)
    % Eqn on p4 of Schafer and Strimmer 2005
    [~, p] = size(x);
    sx = makeMeanZero(x); sx = makeStdOne(sx); % convert S to R
    [r, ~] = varcov(sx);
    lambda = input_lambda;
    lambda = min(lambda, 1); lambda = max(lambda, 0);
    Rhat = (1-lambda)*r;
    Rhat(logical(eye(p))) = 1;
end

function [S, VS] = varcov(x)
    % s(i,j) = cov X(i,j)
    % vs(i,j) = est var s(i,j)
    [n,p] = size(x);
    xc = makeMeanZero(x); 
    S = cov(xc);
    XC1 = repmat(reshape(xc', [p 1 n]), [1 p 1]); % size p*p*n !
    XC2 = repmat(reshape(xc', [1 p n]),  [p 1 1]); % size p*p*n !
    VS = var(XC1 .* XC2, 0,  3) * n/((n-1)^2);
end

function xc = makeMeanZero(x)
% make column means zero
    [n,~] = size(x);
    m = mean(x);
    xc = x - ones(n, 1)*m; 
end

function xc = makeStdOne(x)
    % make column  variances one
    [n,~] = size(x);
    sd = ones(n, 1)*std(x);
    xc = x ./ sd; 
end

function display_graphs(density,char_path,efficiency,mean_cluster_coeff,x_label,title_figure)
    FA_threshold = [0.1:0.1:0.8];

    figure
    hold on 
    plot(FA_threshold,char_path)
    hold on 
    plot(FA_threshold,efficiency)
    xlabel(x_label)
    ylabel('Value')
    legend('mean shortest path','efficiency')
    title(title_figure)

    figure
    plot(FA_threshold,density)
    hold on 
    plot(FA_threshold,mean_cluster_coeff)
    xlabel(x_label)
    ylabel('Value')
    legend('edge density','mean clustering coefficient')
    title(title_figure)
end

function [density,char_path,efficiency,mean_cluster_coeff] = graph_metrics(input)
    % this function return density, number of vertices, and number of edges in 3
    [density,~,~] = density_und(input); 

    % mean shortest path
    % distance matrices 
    distance_mat = distance_bin(input); 
    % this function return network characteristic path length and network global efficiency
    % Characteristic path length is defined here as the mean shortest path length between all pairs of nodes
    % About missing edges and disconnected, it will be infinite value, it
    % would be zero as third argument setting
    [char_path, ~]  = charpath(distance_mat,0,0); 

    % efficiency
    efficiency = 1 / efficiency_bin(input); 
    
    % mean clustering coefficient
    mean_cluster_coeff = mean(clustering_coef_bu(input)); 
end


function show_difference_graphs(input1,input2,y_label)
    threshold = [0.1:0.1:0.8];

    figure
    hold on 
    plot(threshold,input1)
    hold on 
    plot(threshold,input2)
    xlabel('Lambda')
    ylabel(y_label)
    legend('retain negative correlations','discard negative correlations')
    title('effects of retaining or discarding negative correlations')
end