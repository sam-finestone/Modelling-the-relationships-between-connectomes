%% Load data
% Add helper functions to path
addpath('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/Toolboxes/2016_01_16_BCT/');
addpath('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/Toolboxes/NIfTI_20140122/');
addpath('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/Toolboxes/covshrink-kpm/');

% Load csv files
s = zeros(68, 68, 19);
f = zeros(68, 68, 19);
for i = 1:19
    s(:,:,i) = csvread(strcat('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/connectomes-data/Task2Data/', int2str(31+i), '_WFA_68.csv'), 1, 0);
    f(:,:,i) = csvread(strcat('/Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/connectomes-data/Task2Data/', int2str(31+i), '_rsfMRI_68.csv') , 1, 0);
end

%% Calculate t matrix
% Calucluate indirect structural connectivity matrix
t = zeros(size(s));

% Iterate through every patient p
for p = 1:size(s,3)
    % Iterate through (and calculate) every element of t
    for t_i = 1:size(s, 1)
        for t_j = 1:size(s, 2)
    
            % Find the greatest minimum weight in all available two-step chains
            greatest_min_weight = -realmax;
            for k = 1:size(s, 1)
                min_weight = min(s(t_i, k, p), s(k, t_j, p));
                if (min_weight > greatest_min_weight && min_weight ~= 0)
                    greatest_min_weight = min_weight;
                end
            end
    
            % Set matrix element of the indirect structural connectivity matrix
            if (greatest_min_weight < 0)
                greatest_min_weight = 0;
            end
            t(t_i, t_j, p) = greatest_min_weight;
        end 
    end 
end 

%% Fit Model 1: f = alpha + beta * s

% Cross validation loop
best_SSE1 = realmax;
all_alpha1 = zeros(size(s,1), size(s,2), 19);
all_beta1 = zeros(size(s,1), size(s,2), 19);
all_res1 = zeros(size(s,1), size(s,2), 19);

best_alpha1 = zeros(19, 1);
best_beta1  = zeros(19, 1);

best_aic1 = 0; 
best_bic1 = 0; 
for c = 1:19
    % Split data into training/validation set
    training_set = s(:,:,1:size(s,3)~=c);
    training_set_targets = f(:,:,1:size(f,3)~=c);
    validation_set = s(:,:,c);
    validation_set_targets = f(:,:,c);


    alpha1 = 0;
    beta1 = 0; 
    
    % Iterate over every element of the alpha/beta matrix
    
    % Put variables into matrix form: Y = coefficients * X
    X = [ones(68 * 68 * 18, 1) reshape(training_set(:,:,:), 68 * 68 * 18, 1)];
    Y = reshape(training_set_targets(:,:,:), 68 * 68 * 18, 1);

    % Solve for the alpha and beta coefficients
    coeff = pinv(X' * X) * X' * Y;
    alpha1 = coeff(1);
    beta1 = coeff(2);

    % Calculate sum of square errors
    model1_sum = zeros(68);
    temp1 = (alpha1+beta1.*validation_set(:,:));
    res1 = validation_set_targets(:,:)-temp1;
    model1_sum = model1_sum + res1.^2;
    SSE = sum(sum(model1_sum));
    %disp(SSE);
    % Calculate AIC/BIC
    aic1 = 2*3 + 19.*log(model1_sum.*1/19);
    bic1 = 3*log(19) + 19.*log(model1_sum.*1/19);
    % Save the best model
    if SSE < best_SSE1
        best_SSE1 = SSE;
        best_aic1 = mean(mean(aic1));
        best_bic1  = mean(mean(bic1));
    end
    % store alpha for each patient 
    best_alpha1(c) = alpha1;
    best_beta1(c) = beta1;
    
    all_alpha1(:,:,c) = alpha1;
    all_beta1(:,:,c) = beta1;
    all_res1(:,:,c) = res1;
    

end
fprintf('BEST BIC, Model 1: %f\n',best_bic1);
fprintf('BEST AIC, Model 1: %f\n',best_aic1);

%% Fit Model 2: f = alpha + beta * s + y * s^2
% Cross validation loop
best_SSE2 = realmax;

best_alpha2 = zeros(19,1);
best_beta2 = zeros(19,1);
best_y2 = zeros(19,1);

best_bic2 =  0;
best_aic2 =  0; 
for c = 1:19
    % Split data into training/validation set
    training_set = s(:,:,1:size(s,3)~=c);
    training_set_targets = f(:,:,1:size(f,3)~=c);
    validation_set = s(:,:,c);
    validation_set_targets = f(:,:,c);


    alpha2 = 0; 
    beta2 = 0; 
    y2 = 0; 
    
    % Iterate over every element of the alpha/beta/y matrix

    % Put variables into matrix form: Y = coefficients * X
    s_slice = reshape(training_set(:,:,:), 68 * 68 * 18, 1);
    X = [ones(68 * 68 * 18, 1) s_slice (s_slice.^2)];
    Y = reshape(training_set_targets(:,:,:), 68 * 68 * 18, 1);

    % Solve for the alpha and beta coefficients
    coeff = pinv(X' * X) * X' * Y;
    alpha2 = coeff(1);
    beta2 = coeff(2);
    y2 = coeff(3);
    
    % store alpha for each patient 
    best_alpha2(c) = alpha2;
    best_beta2(c) = beta2;
    best_y2(c) = y2;
    
    model2_sum = zeros(68);
    res2 = validation_set_targets(:,:) - (alpha2 + beta2.*validation_set(:,:) + y2.*validation_set(:,:).^2);
    model2_sum = model2_sum + res2.^2;
    SSE = sum(sum(model2_sum));
    % disp(SSE);
    aic2 = 2*4 + 19.*log(model2_sum./19);
    bic2 = 4*log(19) + 19.*log(model2_sum./19);
    % Save the best model
    if SSE < best_SSE2
        best_SSE2 = SSE;
%         best_alpha2 = alpha2;
%         best_beta2 = beta2;
%         best_y2 = y2;
        best_aic2 = mean(mean(aic2));
        best_bic2  = mean(mean(bic2));
    end


end
fprintf('BEST BIC, Model 2: %f\n',best_bic2);
fprintf('BEST AIC, Model 2: %f\n',best_aic2);


%% Fit Model 3: f = alpha + beta * t
% Cross validation loop
best_SSE3 = realmax;

best_alpha3 = zeros(19, 1);
best_beta3 = zeros(19,1);

best_aic3 = 0;
best_bic3  = 0;
for c = 1:19
    % Split data into training/validation set
    training_set = t(:,:,1:size(t,3)~=c);
    training_set_targets = f(:,:,1:size(f,3)~=c);
    validation_set = t(:,:,c);
    validation_set_targets = f(:,:,c);

    alpha3 = 0; % 68 x 68
    beta3 = 0; % 68 x 68
    

    % Put variables into matrix form: Y = coefficients * X
    X = [ones(68 * 68 * 18, 1) reshape(training_set(:,:,:), 68 * 68 * 18, 1)];
    Y = reshape(training_set_targets(:,:,:), 68 * 68 * 18, 1);

    % Solve for the alpha and beta coefficients
    coeff = inv(X' * X) * X' * Y;
    alpha3 = coeff(1);
    beta3 = coeff(2);

    
    % There seems to be some NaN values in the alpha3 and beta3 variables,
    % so these need to be set to zero
    alpha3(isnan(alpha3)) = 0;
    beta3(isnan(beta3)) = 0;
    
    best_alpha3(c) = alpha3;
    best_beta3(c) = beta3;
    
    % Calculate sum of squares
    model3_sum = zeros(68);
    res3 = validation_set_targets(:,:) - (alpha3 + beta3.*validation_set(:,:));
    model3_sum = model3_sum + res3.^2;
    model3_sum(isinf(model3_sum)) = 0; % Set inifity values to zeros 
    SSE = sum(sum(model3_sum));
    %disp(SSE);
    aic3 = 2*3 + 19.* log(model3_sum./19);
    bic3 = 3*log(19) + 19.*log(model3_sum./19);
    % Save the best model
    if SSE < best_SSE3
        best_SSE3 = SSE;
%         best_alpha3 = alpha3;
%         best_beta3 = beta3;
        best_aic3 = mean(mean(aic3));
        best_bic3  = mean(mean(bic3));
    end
    
end
fprintf('BEST BIC, Model 3: %f\n',best_bic3);
fprintf('BEST AIC, Model 3: %f\n',best_aic3);
%% Fit Model 4: f = alpha + beta * t + y * t^2
% Cross validation loop
best_SSE4 = realmax;

% best_alpha4 = zeros(size(s,1), size(s,2));
best_alpha4 = zeros(19,1);
best_beta4 = zeros(19,1);
best_y4 = zeros(19,1);

best_aic4 = 0; 
best_bic4 = 0; 
for c = 1:19
    % Split data into training/validation set
    training_set = t(:,:,1:size(t,3)~=c);
    training_set_targets = f(:,:,1:size(f,3)~=c);
    validation_set = t(:,:,c);
    validation_set_targets = f(:,:,c);

    alpha4 = 0; % 68 x 68
    beta4 = 0; % 68 x 68
    y4 = 0; % 68 x 68
    
    % Iterate over every element of the alpha/beta/y matrix

    % Put variables into matrix form: Y = coefficients * X
    t_slice = reshape(training_set(:,:,:), 68 * 68 * 18, 1);
    X = [ones(68 * 68 * 18, 1) t_slice (t_slice.^2)];
    Y = reshape(training_set_targets(:,:,:), 68 * 68 * 18, 1);
    X(isinf(X)) = 0.00000001; % Set inifity values to zeros 
    X(isnan(X)) = 0.00000001;
    % Solve for the alpha and beta coefficients
%     coeff = pinv(X' * X) * X' * Y;
    coeff = X\Y;
    alpha4 = coeff(1);
    beta4 = coeff(2);
    y4 = coeff(3);


    % There seems to be some NaN values in the alpha3 and beta3 variables,
    % so these need to be set to zero
    alpha4(isnan(alpha4)) = 0;
    beta4(isnan(beta4)) = 0;
    y4(isnan(y4)) = 0;
    
    best_alpha4(c) = alpha4;
    best_beta4(c) = beta4;
    best_y4(c) = y4;
    
    % Calculate sum of square errors
    model4_sum = zeros(68);
    res4 = validation_set_targets(:,:) - (alpha4 + beta4.*validation_set(:,:) + y4.*validation_set(:,:).^2);
    model4_sum = model4_sum + res4.^2;
    model4_sum(isnan(model4_sum)) = 0; % Set NaN values to zeros again 
    model4_sum(isinf(model4_sum)) = 0; % Set inifity values to zeros 
    SSE = sum(sum(model4_sum));
    %disp(SSE);
    
    % AIC/BIC calculation
    aic4 = 2*4 + 19.* log(model4_sum./19);
    bic4 = 4*log(19) + 19.*log(model4_sum./19);

    % Save the best model
    if SSE < best_SSE4
        best_SSE4 = SSE;
%         best_alpha4 = alpha4;
%         best_beta4 = beta4;
%         best_y4 = y4;
        best_aic4 = mean(mean(aic4));
        best_bic4 = mean(mean(bic4));
    end
   
end
fprintf('BEST BIC, Model 4: %f\n',best_bic4);
fprintf('BEST AIC, Model 4: %f\n',best_aic4);

%% Fit Model 5: f = alpha + beta * s + y * t
% Cross validation loop
best_SSE5 = realmax;

best_alpha5 = zeros(19,1);
best_beta5 = zeros(19,1);
best_y5 = zeros(19,1);

best_aic5 = 0;
best_bic5 = 0;
for c = 1:19
    % Split data into training/validation set
    training_set_s = s(:,:,1:size(t,3)~=c);
    training_set_t = t(:,:,1:size(t,3)~=c);
    training_set_targets = f(:,:,1:size(f,3)~=c);
    validation_set_s = s(:,:,c);
    validation_set_t = t(:,:,c);
    validation_set_targets = f(:,:,c);

    alpha5 = 0;
    beta5 = 0;
    y5 = 0;

    % Put variables into matrix form: Y = coefficients * X
    X = [ones(68 * 68 * 18, 1) reshape(training_set_s(:,:,:), 68 * 68 * 18, 1) reshape(training_set_t(:,:,:), 68 * 68 * 18, 1)];
    X(X == 0) = 0.000000001;
    Y = reshape(training_set_targets(:,:,:), 68 * 68 * 18, 1);

    % Solve for the alpha and beta coefficients
    coeff = inv(X' * X) * X' * Y;
    alpha5 = coeff(1);
    beta5 = coeff(2);
    y5 = coeff(3);
    
    % There seems to be some NaN values in the alpha3 and beta3 variables,
    % so these need to be set to zero
    alpha5(isnan(alpha5)) = 0;
    beta5(isnan(beta5)) = 0;
    y5(isnan(y5)) = 0;
    
    best_alpha5(c) = alpha5;
    best_beta5(c) = beta5;
    best_y5(c) = y5;
    
    % Calculate sum of square errors
    model5_sum = zeros(68);
    res5 = validation_set_targets(:,:) - (alpha5 + beta5.*validation_set_s(:,:) + y5.*validation_set_t(:,:));
    model5_sum = model5_sum + res5.^2;
    model5_sum(isnan(model5_sum)) = 0; % Set NaN values to zeros again 
    model5_sum(isinf(model5_sum)) = 0; % Set inifity values to zeros 
    SSE = sum(sum(model5_sum));
    %disp(SSE);
    
    % AIC/BIC
    aic5 = 2*4 + 19.* log(model5_sum./19);
    bic5 = 4*log(19) + 19.*log(model5_sum./19);
    
    % Save the best model
    if SSE < best_SSE5
        best_SSE5 = SSE;
%         best_alpha5 = alpha5;
%         best_beta5 = beta5;
        best_aic5 = mean(mean(aic5));
        best_bic5 = mean(mean(bic5));
    end

end
fprintf('BEST BIC, Model 5: %f\n',best_bic5);
fprintf('BEST AIC, Model 5: %f\n',best_aic5);

disp("Error:")
disp("  Model 1 :" + best_SSE1)
disp("  Model 2 :" + best_SSE2)
disp("  Model 3 :" + best_SSE3)
disp("  Model 4 :" + best_SSE4)
disp("  Model 5 :" + best_SSE5)

disp("AIC:")
disp("  Model 1 :" + best_aic1)
disp("  Model 2 :" + best_aic2)
disp("  Model 3 :" + best_aic3)
disp("  Model 4 :" + best_aic4)
disp("  Model 5 :" + best_aic5)

disp("BIC:")
disp("  Model 1 :" + best_bic1)
disp("  Model 2 :" + best_bic2)
disp("  Model 3 :" + best_bic3)
disp("  Model 4 :" + best_bic4)
disp("  Model 5 :" + best_bic5)

% Plot variance of parameters for the best model (Model 1)
% alpha
mean_alpha1 = mean(best_alpha1);
mean_alpha2 = mean(best_alpha2);
mean_alpha3 = mean(best_alpha3);
mean_alpha4 = mean(best_alpha4);
mean_alpha5 = mean(best_alpha5);

var_alpha1 = var(best_alpha1);
var_alpha2 = var(best_alpha2);
var_alpha3 = var(best_alpha3);
var_alpha4 = var(best_alpha4);
var_alpha5 = var(best_alpha5);

disp("Means alpha:")
disp("  Model 1 :" + mean_alpha1)
disp("  Model 2 :" + mean_alpha2)
disp("  Model 3 :" + mean_alpha3)
disp("  Model 4 :" + mean_alpha4)
disp("  Model 5 :" + mean_alpha5)

disp("Variance alpha:");
disp("  Model 1 :" + var_alpha1)
disp("  Model 2 :" + var_alpha2)
disp("  Model 3 :" + var_alpha3)
disp("  Model 4 :" + var_alpha4)
disp("  Model 5 :" + var_alpha5)

% beta parameter 
mean_beta1 = mean(best_beta1);
mean_beta2 = mean(best_beta2);
mean_beta3 = mean(best_beta3);
mean_beta4 = mean(best_beta4);
mean_beta5 = mean(best_beta5);

var_beta1 = var(best_beta1);
var_beta2 = var(best_beta2);
var_beta3 = var(best_beta3);
var_beta4 = var(best_beta4);
var_beta5 = var(best_beta5);

disp("Means beta:")
disp("  Model 1 :" + mean_beta1)
disp("  Model 2 :" + mean_beta2)
disp("  Model 3 :" + mean_beta3)
disp("  Model 4 :" + mean_beta4)
disp("  Model 5 :" + mean_beta5)

disp("Variance beta:");
disp("  Model 1 :" + var_beta1)
disp("  Model 2 :" + var_beta2)
disp("  Model 3 :" + var_beta3)
disp("  Model 4 :" + var_beta4)
disp("  Model 5 :" + var_beta5)

% gamma paramter 

mean_gamma2 = mean(best_y2);
mean_gamma4 = mean(best_y4);
mean_gamma5 = mean(best_y5);

var_gamma2 = var(best_y2);
var_gamma4 = var(best_y4);
var_gamma5 = var(best_y5);

disp("Means gamma:")
disp("  Model 2 :" + mean_gamma2)
disp("  Model 4 :" + mean_gamma4)
disp("  Model 5 :" + mean_gamma5)

disp("Variance gamma:");
disp("  Model 2 :" + var_gamma2)
disp("  Model 4 :" + var_gamma4)
disp("  Model 5 :" + var_gamma5)

% subplot(2,3,1);
% imshow(var(all_alpha1,0,3) / max(max(var(all_alpha1,0,3))));
% title("Alpha Variance");
% subplot(2,3,2);
% imshow(var(all_beta1,0,3) / max(max(var(all_beta1,0,3))));
% title("Beta Variance");
% % Plot the average sum of square error for the best model (Model 1)
% subplot(2,3,3);
% imshow(mean(all_res1,3) / max(max(mean(all_res1,3))));
% title("Performance");

% Plot variance of parameters for the best model (Model 1)
% subplot(2, 3, 4);
% zoom_start_x = 20;
% zoom_end_x = 30;
% zoom_start_y = 10;
% zoom_end_y = 20;
% imshow(var(all_alpha1(zoom_start_x:zoom_end_x, zoom_start_y:zoom_end_y,:),0,3) / max(max(var(all_alpha1(zoom_start_x:zoom_end_x, zoom_start_y:zoom_end_y,:),0,3))));
% title("Zoomed Alpha Variance");
% subplot(2,3,5);
% imshow(var(all_beta1(zoom_start_x:zoom_end_x, zoom_start_y:zoom_end_y,:),0,3) / max(max(var(all_beta1(zoom_start_x:zoom_end_x, zoom_start_y:zoom_end_y,:),0,3))));
% title("Zoomed Beta Variance");
% % Plot the average sum of square error for the best model (Model 1)
% subplot(2,3,6);
% imshow(mean(all_res1(zoom_start_x:zoom_end_x, zoom_start_y:zoom_end_y,:),3) / max(max(mean(all_res1(zoom_start_x:zoom_end_x, zoom_start_y:zoom_end_y,:),3))));
% title("Zoomed Performance");

function m = MeanExcludingInf(matrix)
    
    sum = 0;
    count = 0;
    for i=1:size(matrix,1)
        for j=1:size(matrix,2)
            if isinf(matrix(i,j)) == 0
                sum = sum + matrix(i,j);d
                count = count + 1;
            end
        end
    end
    m = sum / count;
end
