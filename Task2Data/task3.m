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


%% Task 1 
structural_vertex_density = zeros(19, 68);
functional_vertex_density = zeros(19, 68);

% Estimate the structural connectivity density for each vertex
for i = 1:19
    % load the structural connectomes 
    temp = sum(s(:,:,i), 1); 
    structural_vertex_density(i,:) = sum(s(:,:,i), 1); % same as rows
    functional_vertex_density(i,:) = sum(f(:,:,i), 1); % same as rows
end

% model functional connectivity density based on structural connectivity 
% density, independently for each vertex. on a single fold
% => f = alpha + beta * s
alphas = zeros(68,1);
betas = zeros(68,1);
gammas = zeros(68,1);
j = 1;
SSE_ = zeros(68, 1);

training_set = structural_vertex_density(1:size(structural_vertex_density,1)~=j, :);
training_set_targets = functional_vertex_density(1:size(functional_vertex_density,1)~=j, :);

validation_set = structural_vertex_density(j, :);
validation_set_targets = functional_vertex_density(j, :);

% model functional connectivity density based on structural connectivity
% density, independently for each vertex.
for i=1:68
    feature = reshape(training_set(:,i),  18*1, 1);
%     X = [ones(1 * 18, 1) feature (feature.^2)];
    X = [ones(1 * 18, 1) feature];
    Y = reshape(training_set_targets(:, i), 18*1, 1);
    % Solve for the alpha and beta coefficients
    coeff = pinv(X' * X) * X' * Y;
    alpha = coeff(1);
    beta = coeff(2);
%     gamma = coeff(3);
    % Calculate sum of square errors
%     model1_sum = zeros(68);
    res1 = validation_set_targets(:,i)-(alpha+beta.*validation_set(:,i));
%     res1 = validation_set_targets(:,i) - (alpha + beta.*validation_set(:,i) + gamma.*validation_set(:,i).^2);
    model1_sum = res1.^2;
    SSE = sum(model1_sum);
    SSE_(i) = SSE;
end

% mse
fprintf('---- mean squared error in holdout set ---- \n');
fprintf('MSE, OLS: %f\n',SSE);

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
            if (greatest_min_weight < 0)
                greatest_min_weight = 0;
            end
            % Set matrix element of the indirect structural connectivity matrix
            t(t_i, t_j, p) = greatest_min_weight;
        end 
    end 
end 
%% Fit Model 1: f = alpha + beta * s

% Cross validation loop
SSEs =  zeros(19, 1);
best_SSE1 = realmax;

best_alpha1 = zeros(19, 1);
best_beta1  = zeros(19, 1);

best_aic1 = 0; 
best_bic1 = 0; 
for c = 1:19
    % Split data into training/validation set
    training_set = s(:,:,c);
    training_set_targets = f(:,:,c);
    n = [1:19];
    n(n==c) = [];
    index_rand = randi(18);
    validation_set = s(:,:,index_rand);
    validation_set_targets = f(:,:,index_rand);


    alpha1 = 0;
    beta1 = 0; 
    
    % Iterate over every element of the alpha/beta matrix
    
    % Put variables into matrix form: Y = coefficients * X
    X = [ones(68 * 68 * 1, 1) reshape(training_set(:,:,:), 68 * 68 * 1, 1)];
    Y = reshape(training_set_targets(:,:,:), 68 * 68 * 1, 1);

    % Solve for the alpha and beta coefficients
    coeff = pinv(X' * X) * X' * Y;
    alpha1 = coeff(1);
    beta1 = coeff(2);
    % store alpha for each patient 
    best_alpha1(c) = alpha1;
    best_beta1(c) = beta1;
    
    % Calculate sum of square errors
    model1_sum = zeros(68);
    res1 = validation_set_targets(:,:)-(alpha1+beta1.*validation_set(:,:));
    model1_sum = model1_sum + res1.^2;
    SSE = sum(sum(model1_sum));
    SSEs(c) = SSE;
    %disp(SSE);
    % Calculate AIC/BIC
    aic1 = 2*3 + 19.*log(model1_sum.*1/19);
    bic1 = 3*log(19) + 19.*log(model1_sum.*1/19);
    % store alpha for each patient 
    % Save the best model
    if SSE < best_SSE1
        best_SSE1 = SSE;
        best_aic1 = mean(mean(aic1));
        best_bic1  = mean(mean(bic1));
    end
    
end

disp("SSE:")
% eval model
disp("  Model 1 :" + SSEs);
fprintf('MSE, Model 1: %f\n',mean(SSEs));
fprintf('BEST BIC, Model 1: %f\n',best_bic1);
fprintf('BEST AIC, Model 1: %f\n',best_aic1);
% eval param 
fprintf('MEAN ALPHAS, Model 1: %f\n',mean(best_alpha1));
fprintf('MEAN VARIANCE, Model 1: %f\n',mean(best_beta1));
fprintf('VARIANCE ALPHAS, Model 1: %f\n',var(best_alpha1));
fprintf('VARIANCE BETAS, Model 1: %f\n',var(best_beta1));
%% Fit Model 2: f = alpha + beta * s + y * s^2
% Cross validation loop
SSE2s = zeros(19, 1);
best_SSE2 = realmax;

best_alpha2 = zeros(19,1);
best_beta2 = zeros(19,1);
best_y2 = zeros(19,1);

best_bic2 =  0;
best_aic2 =  0; 
for c = 1:19
    % Split data into training/validation set
    training_set = s(:,:,c);
    training_set_targets = f(:,:,c);
    n = [1:19];
    n(n==c) = [];
    index_rand = randi(18);
    validation_set = s(:,:,index_rand);
    validation_set_targets = f(:,:,index_rand);

    alpha2 = 0; 
    beta2 = 0; 
    y2 = 0; 
    
    % Iterate over every element of the alpha/beta/y matrix

    % Put variables into matrix form: Y = coefficients * X
    s_slice = reshape(training_set(:,:,:), 68 * 68 * 1, 1);
    X = [ones(68 * 68 * 1, 1) s_slice (s_slice.^2)];
    Y = reshape(training_set_targets(:,:,:), 68 * 68 * 1, 1);

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
    SSE2s(c) = SSE;
    % disp(SSE);
    aic2 = 2*4 + 19.*log(model2_sum./19);
    bic2 = 4*log(19) + 19.*log(model2_sum./19);
   if SSE < best_SSE2
        best_SSE2 = SSE;
        best_aic2 = mean(mean(aic2));
        best_bic2  = mean(mean(bic2));
    end
    
end

disp("SSE:")
disp("Model 2 :" + SSE2s);
fprintf('MSE, Model 2: %f\n',mean(SSE2s));
fprintf('BEST BIC, Model 2: %f\n',best_bic2);
fprintf('BEST AIC, Model 2: %f\n',best_aic2);

% eval param 
fprintf('MEAN ALPHAS, Model 2: %f\n',mean(best_alpha2));
fprintf('MEAN BETAS, Model 2: %f\n',mean(best_beta2));
fprintf('MEAN GAMMAS, Model 2: %f\n',mean(best_y2));
fprintf('VARIANCE ALPHAS, Model 2: %f\n',var(best_alpha2));
fprintf('VARIANCE BETAS, Model 2: %f\n',var(best_beta2));
fprintf('VARIANCE GAMMAS, Model 2: %f\n',var(best_y2));

%% Fit Model 3: f = alpha + beta * t
% Cross validation loop
SSE3s = zeros(19,1);
best_SSE3 = realmax;

best_alpha3 = zeros(19, 1);
best_beta3 = zeros(19,1);

best_aic3 = 0;
best_bic3  = 0;
for c = 1:19
    % Split data into training/validation set
    training_set = s(:,:,c);
    training_set_targets = f(:,:,c);
    n = [1:19];
    n(n==c) = [];
    index_rand = randi(18);
    validation_set = s(:,:,index_rand);
    validation_set_targets = f(:,:,index_rand);

    alpha3 = 0; % 68 x 68
    beta3 = 0; % 68 x 68
    

    % Put variables into matrix form: Y = coefficients * X
    X = [ones(68 * 68 * 1, 1) reshape(training_set(:,:,:), 68 * 68 * 1, 1)];
    Y = reshape(training_set_targets(:,:,:), 68 * 68 * 1, 1);

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
    SSE3s(c) = SSE;
    %disp(SSE);
    aic3 = 2*3 + 19.* log(model3_sum./19);
    bic3 = 3*log(19) + 19.*log(model3_sum./19);
    % Save the best model
    if SSE < best_SSE3
        best_SSE3 = SSE;
        best_aic3 = mean(mean(aic3));
        best_bic3  = mean(mean(bic3));
    end
    
end

disp("SSE:")
disp("  Model 3:" + SSE3s);
fprintf('MSE, Model 3: %f\n',mean(SSE3s));
fprintf('BEST BIC, Model 3: %f\n',best_bic3);
fprintf('BEST AIC, Model 3: %f\n',best_aic3);
% eval param 
fprintf('MEAN ALPHAS, Model 3: %f\n',mean(best_alpha3));
fprintf('MEAN VARIANCE, Model 3: %f\n',mean(best_beta3));
fprintf('VARIANCE ALPHAS, Model 3: %f\n',var(best_alpha3));
fprintf('VARIANCE BETAS, Model 3: %f\n',var(best_beta3));
%% Fit Model 4: f = alpha + beta * t + y * t^2
% Cross validation loop
SSE4s = zeros(19,1);
best_SSE4 = realmax;

% best_alpha4 = zeros(size(s,1), size(s,2));
best_alpha4 = zeros(19,1);
best_beta4 = zeros(19,1);
best_y4 = zeros(19,1);

best_aic4 = 0; 
best_bic4 = 0; 
for c = 1:19
    % Split data into training/validation set
    training_set = s(:,:,c);
    training_set_targets = f(:,:,c);
    n = [1:19];
    n(n==c) = [];
    index_rand = randi(18);
    validation_set = s(:,:,index_rand);
    validation_set_targets = f(:,:,index_rand);

    alpha4 = 0; % 68 x 68
    beta4 = 0; % 68 x 68
    y4 = 0; % 68 x 68
    
    % Iterate over every element of the alpha/beta/y matrix

    % Put variables into matrix form: Y = coefficients * X
    t_slice = reshape(training_set(:,:,:), 68 * 68 * 1, 1);
    X = [ones(68 * 68 * 1, 1) t_slice (t_slice.^2)];
    Y = reshape(training_set_targets(:,:,:), 68 * 68 * 1, 1);
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
    SSE4s(c) = SSE;
    %disp(SSE);
    
    % AIC/BIC calculation
    aic4 = 2*4 + 19.* log(model4_sum./19);
    bic4 = 4*log(19) + 19.*log(model4_sum./19);

    % Save the best model
    if SSE < best_SSE4
        best_SSE4 = SSE;
        best_aic4 = mean(mean(aic4));
        best_bic4 = mean(mean(bic4));
    end
end
disp("SSE:")
disp("Model 4:" + SSE4s);
fprintf('MSE, Model 4: %f\n',mean(SSE4s));
fprintf('BEST BIC, Model 4: %f\n',best_bic4);
fprintf('BEST AIC, Model 4: %f\n',best_aic4);

% eval param 
fprintf('MEAN ALPHAS, Model 4: %f\n',mean(best_alpha4));
fprintf('MEAN BETAS, Model 4: %f\n',mean(best_beta4));
fprintf('MEAN GAMMAS, Model 4: %f\n',mean(best_y4));
fprintf('VARIANCE ALPHAS, Model 4: %f\n',var(best_alpha4));
fprintf('VARIANCE BETAS, Model 4: %f\n',var(best_beta4));
fprintf('VARIANCE GAMMAS, Model 4: %f\n',var(best_y4));

% plot the relationship
% figure; 
% for i=1:10   
%     position = randi([1,68],2,1);
%     fval = [];
%     sval = [];
%     for k = 1:19
%         fval = [fval,f(position(1),position(2),k)];
%         sval = [sval,s(position(1),position(2),k)];
%     end
%     
%     subplot(2,5,i);
%     x = linspace(0,0.7);
%     a = best_alpha4(position(1));
%     b = best_beta4(position(1));
%     c = best_y4(position(1));
%     plot(x,a+b.*x+c.*x.^2);
%     hold on;
%     
%     scatter(sval,fval);
%     title(sprintf('(%.0f , %.0f)', position(1),position(2)));
% end
%% Fit Model 5: f = alpha + beta * s + y * t
% Cross validation loop
SSE5s = zeros(19,1);
best_SSE5 = realmax;

best_alpha5 = zeros(19,1);
best_beta5 = zeros(19,1);
best_y5 = zeros(19,1);

best_aic5 = 0;
best_bic5 = 0;
for c = 1:19
    % Split data into training/validation set
    training_set_s = s(:,:,c);
    training_set_targets = f(:,:,c);
    n = [1:19];
    n(n==c) = [];
    index_rand = randi(18);
    validation_set_s = s(:,:,index_rand);
    validation_set_targets = f(:,:,index_rand);
    training_set_t = t(:,:,c);
    validation_set_t = t(:,:,index_rand);

    alpha5 = 0;
    beta5 = 0;
    y5 = 0;

    % Put variables into matrix form: Y = coefficients * X
    X = [ones(68 * 68 * 1, 1) reshape(training_set_s(:,:,:), 68 * 68 * 1, 1) reshape(training_set_t(:,:,:), 68 * 68 * 1, 1)];
    X(X == 0) = 0.000000001;
    Y = reshape(training_set_targets(:,:,:), 68 * 68 * 1, 1);

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
    SSE5s(c) = SSE;
    %disp(SSE);
    
    % AIC/BIC
    aic5 = 2*4 + 19.* log(model5_sum./19);
    bic5 = 4*log(19) + 19.*log(model5_sum./19);
    
    % Save the best model
    if SSE < best_SSE5
        best_SSE5 = SSE;
        best_aic5 = mean(mean(aic5));
        best_bic5 = mean(mean(bic5));
    end
end

disp("SSE:");
disp("Model 5:" + SSE5s);
fprintf('MSE, Model 5: %f\n',mean(SSE5s)); 
fprintf('BEST BIC, Model 5: %f\n',best_bic5);
fprintf('BEST AIC, Model 5: %f\n',best_aic5);
% eval param 
fprintf('MEAN ALPHAS, Model 5: %f\n',mean(best_alpha5));
fprintf('MEAN BETAS, Model 5: %f\n',mean(best_beta5));
fprintf('MEAN GAMMAS, Model 5: %f\n',mean(best_y5));
fprintf('VARIANCE ALPHAS, Model 5: %f\n',var(best_alpha5));
fprintf('VARIANCE BETAS, Model 5: %f\n',var(best_beta5));
fprintf('VARIANCE GAMMAS, Model 5: %f\n',var(best_y5));
%% Model 6: f = alpha + beta * s + gamma * t.^2
% Cross validation loop
SSE6s = zeros(19,1);
best_SSE6 = realmax;

best_alpha6 = zeros(19,1);
best_beta6 = zeros(19,1);
best_y6 = zeros(19,1);

best_aic6 = 0;
best_bic6 = 0;
for c = 1:19
    % Split data into training/validation set
    training_set_s = s(:,:,c);
    training_set_targets = f(:,:,c);
    n = [1:19];
    n(n==c) = [];
    index_rand = randi(18);
    validation_set_s = s(:,:,index_rand);
    validation_set_targets = f(:,:,index_rand);
    training_set_t = t(:,:,c);
    validation_set_t = t(:,:,index_rand);

    alpha6 = 0;
    beta6 = 0;
    y6 = 0;

    % Put variables into matrix form: Y = coefficients * X
    X = [ones(68 * 68 * 1, 1) reshape(training_set_s(:,:,:), 68 * 68 * 1, 1) reshape(training_set_t(:,:,:), 68 * 68 * 1, 1)];
    X(X == 0) = 0.000000001;
    Y = reshape(training_set_targets(:,:,:), 68 * 68 * 1, 1);

    % Solve for the alpha and beta coefficients
    coeff = inv(X' * X) * X' * Y;
    alpha6 = coeff(1);
    beta6 = coeff(2);
    y6 = coeff(3);
    
    % There seems to be some NaN values in the alpha3 and beta3 variables,
    % so these need to be set to zero
    alpha6(isnan(alpha6)) = 0;
    beta6(isnan(beta6)) = 0;
    y6(isnan(y6)) = 0;
    
    best_alpha6(c) = alpha6;
    best_beta6(c) = beta6;
    best_y6(c) = y6;
    
    % Calculate sum of square errors
    model6_sum = zeros(68);
    res6 = validation_set_targets(:,:) - (alpha6 + beta6.*validation_set_t(:,:) + y6.*validation_set_s(:,:).^2);
    model6_sum = model6_sum + res6.^2;
    model6_sum(isnan(model6_sum)) = 0; % Set NaN values to zeros again 
    model6_sum(isinf(model6_sum)) = 0; % Set inifity values to zeros 
    SSE = sum(sum(model6_sum));
    SSE6s(c) = SSE;
    %disp(SSE);
    
    % AIC/BIC
    aic6 = 2*4 + 19.* log(model6_sum./19);
    bic6 = 4*log(19) + 19.*log(model6_sum./19);
    
    % Save the best model
    if SSE < best_SSE6
        best_SSE6 = SSE;
        best_aic6 = mean(mean(aic6));
        best_bic6 = mean(mean(bic6));
    end
    
end

disp("SSE:");
disp("Model 6:" + SSE6s);
% model param 
fprintf('MSE, Model 6: %f\n', mean(SSE6s));
fprintf('BEST BIC, Model 6: %f\n', best_bic6);
fprintf('BEST AIC, Model 6: %f\n', best_aic6);
% eval param 
fprintf('MEAN ALPHAS, Model 6: %f\n',mean(best_alpha6));
fprintf('MEAN BETAS, Model 6: %f\n',mean(best_beta6));
fprintf('MEAN GAMMAS, Model 6: %f\n',mean(best_y6));
fprintf('VARIANCE ALPHAS, Model 6: %f\n',var(best_alpha6));
fprintf('VARIANCE BETAS, Model 6: %f\n',var(best_beta6));
fprintf('VARIANCE GAMMAS, Model 6: %f\n',var(best_y6));

%% Use sparse linear models based on LASSO to 
% relate each functional connection to a subset
% of structural connections f = alpha + beta*s

% Cross validation loop
alpha_lasso_all = zeros(size(s,1), size(s,2), 19);
beta_lasso_all = zeros(size(s,1), size(s,2), 19);
all_res1 = zeros(size(s,1), size(s,2), 19);
best_SSE1 = realmax;

best_alpha_values = zeros(19, 1);
best_beta1  = zeros(19, 1);

best_aic1 = 0; 
best_bic1 = 0; 
avg_SSE_lassso = zeros(19, 1);
k = 3;
indices = crossvalind('Kfold',19,k);
total_sse_lasso = zeros(68, 68, k);

% Do 
% c = randi(19); 
for c = 1:3
%     Split data into training/validation set
    index_train = find(indices~=c);
    index_test = find(indices==c);
    train_size = length(index_train);
    training_set = s(:,:,index_train);
    training_set_targets = f(:,:,index_train);
    validation_set = s(:,:,index_test);
    validation_set_targets = f(:,:,index_test);
% training_set = s(:,:,1:19 ~= c);
% training_set_targets = f(:,:,1:19 ~= c);
% validation_set = s(:,:,c);
% validation_set_targets = f(:,:,c);

    alpha_lasso = zeros(size(s,1),size(s,2)); 
    beta_lasso = zeros(size(s,1),size(s,2)); 

    % Iterate over every element of the alpha/beta matrix
    for i = 1:size(s, 1)
        for j = 1:size(s, 2)
            % Put variables into matrix form: Y = coefficients * X
    %             X = [ones(18, 1) reshape(training_set, 68*68, 18)];
            X = reshape(training_set, 68*68, train_size)';
            z = training_set_targets(i,j,:);
            Y = reshape(training_set_targets(i,j,:), 1, train_size)';

            [Blasso, stats] = lasso(X, Y);
%             lambda = [1:length(index_test)] + 80;
            [~, total_lambda] = size(Blasso);
            lambda = floor(total_lambda*0.7);
            beta_lasso = Blasso(:, lambda)';   
            target = reshape(validation_set_targets(i,j, :), 1, length(index_test))';
            temp = reshape(validation_set, 68*68, length(index_test));
            predicted = beta_lasso*temp;
            
            residual = target - predicted';
            sse = residual.^2;
            total_sse_lasso(i,j, c) = mean(sse);
        end
    end
end 

average_see_lasso = mean(mean(total_sse_lasso));
fprintf('MSE, LASSO: %f\n',average_see_lasso);

csvwrite('output_lasso_avg.csv',average_see_lasso);
csvwrite('output_lasso_tot.csv',total_sse_lasso);

%%
% Split data into training/validation set
%     training_set = s(:,:,1:size(s,3)~=c);
%     training_set_targets = f(:,:,1:size(f,3)~=c);
%     validation_set = s(:,:,c);
%     validation_set_targets = f(:,:,c);
% 
% 
%     alpha_lasso = zeros(size(s,1),size(s,2)); % 68 x 68
%     beta_lasso = zeros(size(s,1),size(s,2)); % 68 x 68
%     
%     % Iterate over every element of the alpha/beta matrix
%     for i = 1:size(s, 1)
%         for j = 1:size(s, 2)
%             % Put variables into matrix form: Y = coefficients * X
% %             X = [ones(18, 1) reshape(training_set, 68*68, 18)];
%             X = reshape(training_set, 68*68, 18)';
%             Y = reshape(training_set_targets(i,j,:), 1, 18)';
%     
%             % Solve for the alpha and beta coefficients
%             [B, stats] = lasso(X,Y);
% %             idxLambdaMinMSE = FitInfo.IndexMinMSE;
% %             minMSEModelPredictors = FitInfo.PredictorNames(B(:,idxLambdaMinMSE)~=0);
%             coeff_lasso = B(:,1);
%             alpha_lasso(i,j) = coeff_lasso(1);
%             beta_lasso(i,j) = coeff_lasso(2);
%         end
%     end
%     % There seems to be some NaN values in the alpha3 and beta3 variables,
%     % so these need to be set to zero
%     alpha_lasso(isnan(alpha_lasso)) = 0;
%     beta_lasso(isnan(beta_lasso)) = 0;
%     % Calculate sum of square errors for the OLs
%     model1_sum = zeros(68);
%     res_lasso = validation_set_targets(:,:) - (alpha_lasso+beta_lasso.*validation_set(:,:));
%     model1_sum = model1_sum + res_lasso.^2;
%     SSE_lasso = sum(sum(model1_sum));
%     avg_SSE_lassso(c) = SSE_lasso;
%     % Calculate sum of square errors for the LASSO
% %     model1_sum = zeros(68);
% %     res_ols = validation_set_targets(:,:)-(alpha1+beta1.*validation_set(:,:));
% %     model1_sum = model1_sum + res_ols.^2;
% %     SSE_ols = sum(sum(model1_sum));
%     
%     
%     %disp(SSE);
%     % Calculate AIC/BIC
%     aic1 = 2*3 + 19.*log(model1_sum.*1/19);
%     bic1 = 3*log(19) + 19.*log(model1_sum.*1/19);
%     % Save the best model
%     if SSE_lasso < best_SSE1
%         best_SSE1 = SSE_lasso;
%         best_aic1 = mean(mean(aic1));
%         best_bic1  = mean(mean(bic1));
%     end
%     % store alpha for each patient 
% %     best_alpha1(c) = alpha1;
% %     best_beta1(c) = beta1;
%     
%     beta_lasso_all(:,:,c) = beta_lasso;
%     alpha_lasso_all(:,:,c) = alpha_lasso;
%     all_res1(:,:,c) = res_lasso;
%     
% 
% end
% 
% fprintf('MSE, LASSO: %f\n',mean(avg_SSE_lassso));

% Split data into training/validation set
% c = 1;
% training_set = s(:,:,1:size(s,3)~=c);
% training_set_targets = f(:,:,1:size(f,3)~=c);
% validation_set = s(:,:,c);
% validation_set_targets = f(:,:,c);
% 
% 
% alpha1 = 0;
% beta1 = 0; 
% 
% % Iterate over every element of the alpha/beta matrix
% 
% % Put variables into matrix form: Y = coefficients * X
% X = [ones(68 * 68 * 18, 1) reshape(training_set(:,:,:), 68 * 68 * 18, 1)];
% Y = reshape(training_set_targets(:,:,:), 68 * 68 * 18, 1);
% 
% % Solve for the alpha and beta coefficients using ols
% % coeff_ols = pinv(X' * X) * X' * Y;
% coeff_ols = (X' * X) \ X' * Y;
% 
% % Solve for the alpha and beta coeff usisng LASSO
% [coeff_lasso, stats] = lasso(X, Y, 'CV', 18);
% 
% lassoPlot(coeff_lasso, stats, 'PlotType', 'CV');
% Blasso = [stats.Intercept(stats.Index1SE); coeff_lasso(:, stats.Index1SE)];
% 
% % Evaluate on the holdout set 
% X_val = [ones(68 * 68, 1) reshape(validation_set(:,:), 68 * 68, 1)];
% yLASSO = X_val*Blasso;
% yOLS = X_val*coeff_ols;
% 
% % mse
% fprintf('---- mean squared error in holdout set ---- \n');
% fprintf('MSE, OLS: %f\n',mean((validation_set_targets(:,:) - yOLS).^2));
% fprintf('MSE, LASSO: %f\n',mean((validation_set_targets(:,:) - yLASSO).^2));
% 
% alpha1 = coeff_ols(1);
% beta1 = coeff_ols(2);
% alpha_lasso = coeff_lasso(1);
% beta_lasso = coeff_lasso(2);
% 
% % Calculate sum of square errors for the OLs
% model1_sum = zeros(68);
% res_ols = validation_set_targets(:,:) - (alpha1+beta1.*validation_set(:,:));
% model1_sum = model1_sum + res_ols.^2;
% SSE_ols = sum(sum(model1_sum));
% 
% % Calculate sum of square errors for the LASSO
% model1_sum_l = zeros(68);
% res_lasso = validation_set_targets(:,:)-(alpha_lasso+beta_lasso.*validation_set(:,:));
% model1_sum_l = model1_sum_l + res_lasso.^2;
% SSE_lasso = sum(sum(model1_sum_l));
% 
% fprintf('---- sum squared error in holdout set ---- \n');
% fprintf('SSE, OLS: %f\n',SSE_ols);
% fprintf('SSE, LASSO: %f\n',SSE_lasso);
% 
% 
% % Histogram of the estimated coefficients 
% histogram(coeff_lasso, 20)
% hold on;
% histogram(coeff_ols, 20);
% hold off;
% legend('OLS', 'LASSO'); 

% scatter plot 
% xx = linspace(min([coeff_ols: coeff_lasso]), max([coeff_ols: coeff_lasso]));
% plot(coeff_ols, coeff_lasso, 'o', xx, xx, '-');
% legend('OLS', 'lasso');

function m = MeanExcludingInf(matrix)
    
    sum = 0;
    count = 0;
    for i=1:size(matrix,1)
        for j=1:size(matrix,2)
            if isinf(matrix(i,j)) == 0
                sum = sum + matrix(i,j);
                count = count + 1;
            end
        end
    end
    m = sum / count;
end
