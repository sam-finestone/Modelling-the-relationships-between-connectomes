% Estimate the structural connectivity density for each vertex (the sum across rows or
% columns of each structural connectivity matrix). Do the same for the functional connectivity
% matrices, and then model functional connectivity density based on structural connectivity
% density, independently for each vertex. Do you find any strong associations between the
% modalities using this approach?
close all;
clear all;

% Add helper functions to path
addpath('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/Toolboxes/2016_01_16_BC');
addpath('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/Toolboxes/NIfTI_20140122');
addpath('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/Toolboxes/covshrink-kpm');


FA_threshold = [0.1:0.1:0.8];
connectome_file = {'FA.1_graph.csv'; 'FA.2_graph.csv'; 'FA.3_graph.csv';...
            'FA.4_graph.csv'; 'FA.5_graph.csv'; 'FA.6_graph.csv'; ...
            'FA.7_graph.csv'; 'FA.8_graph.csv'};

structural_vertex_density = zeros(8,68);

for i = 1:8
    % load the structural connectomes 
    fa_file_path = ['/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/connectomes-data/Task1Data/connectomes/', connectome_file{i}];
    % skip the first three lines 
    struct_connect{i} = csvread(fa_file_path,3,0);
    structural_vertex_density(i,:) = sum(struct_connect{i}, 1); % same as rows
end


parcellation_func = load_nii('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/connectomes-data/Task1Data/tractor/functional/parcellation.nii.gz');
data_func = load_nii('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/connectomes-data/Task1Data/tractor/functional/data.nii.gz');

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

functional_vertex_density = zeros(8,108);

% resulting matrics at a correlation of 0.1, binarise it
for t = 1:8
    corshrink_new = corshrink_matrix{t};
    % binarise it
    corshrink_new(corshrink_new < 0.1) = 0;
    corshrink_new(corshrink_new >= 0.1) = 1;
    binary_corshrink{t} = corshrink_new;
    functional_vertex_density(t,:) = sum(binary_corshrink{t}, 1); % same as rows
end

% model functional connectivity density based on structural connectivity 
% density, independently for each vertex. 
% => f = alpha + beta * s
for j=1:8
    % Put variables into matrix form: Y = coefficients * X
    X = [ones(68 * 8, 1) reshape(structural_vertex_density, 68 * 8, 1)];
    Y = reshape(functional_vertex_density, 108 * 8, 1);

    % Solve for the alpha and beta coefficients
    coeff = pinv(X' * X) * X' * Y;
    alpha1 = coeff(1);
    beta1 = coeff(2);
end 



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





