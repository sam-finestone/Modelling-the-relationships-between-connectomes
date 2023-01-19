% directory = dir('.');
% check_data = csvread('32_rsfMRI_68.csv' , 1); 
% 
% rsfMRI_data = zeros(19,68,68);
% WFA_data = zeros(19,68,68);
% i = 1;
% j = 1;
% idx = 32;
% for k = 2:length(directory)
%     baseFileName = directory(k).name;
%     if endsWith(baseFileName, strcat('_', 'rsfMRI', '_68.csv'))
% %         disp(baseFileName);
%         data = csvread(baseFileName,1);
%         rsfMRI_data(i,:,:) = data;
%         i = i + 1;
%     end 
%     if endsWith(baseFileName, strcat('_', 'WFA', '_68.csv'))
%         disp(baseFileName);
%         data = csvread(baseFileName,1);
%         WFA_data(j,:,:) = data;
%         j = j + 1;
%     end
%     idx = idx + 1/2;
% end
% 
% 
% % tij = maxk { min { sik, skj } } s.t. sik, skj ≠ 0
% 
% tij = zeros(19,68,68);
% for subject=1:19
%     for i=1:68 
%         for j=1:68
%             % Find the greatest minimum weight in all available two-step chains
%             greatest_min_weight = -realmax;
%             for k = 1:size(s, 1)
%                 min_weight = min(s(t_i, k, p), s(k, t_j, p));
%                 if (min_weight > greatest_min_weight)
%                     greatest_min_weight = min_weight;
%                 end
%             end
%     
%             % Set matrix element of the indirect structural connectivity matrix
%             t(t_i, t_j, p) = greatest_min_weight;
%         end 
%     end 
% end 


% cd /Users/sam/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Masters/UCL/Masters/T2/CMBI/CourseWork/cw3/connectomes-data/Task2Data

% Based on Harvey's code

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
    res1 = validation_set_targets(:,:)-(alpha1+beta1.*validation_set(:,:));
    model1_sum = model1_sum + res1.^2;
    SSE = sum(sum(model1_sum));
    %disp(SSE);

    % Save the best model
    if SSE < best_SSE1
        best_SSE1 = SSE;
    end
    all_alpha1(:,:,c) = alpha1;
    all_beta1(:,:,c) = beta1;
    all_res1(:,:,c) = res1;
    
    % Calculate AIC/BIC
    aic1 = 2*3 + 19.*log(model1_sum.*1/19);
    bic1 = 3*log(19) + 19.*log(model1_sum.*1/19);

end

%% Fit Model 2: f = alpha + beta * s + y * s^2
% Cross validation loop
best_SSE2 = realmax;
best_alpha2 = zeros(size(s,1), size(s,2));
best_beta2 = zeros(size(s,1), size(s,2));
best_y2 = zeros(size(s,1), size(s,2));
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

    
    model2_sum = zeros(68);
    res2 = validation_set_targets(:,:) - (alpha2 + beta2.*validation_set(:,:) + y2.*validation_set(:,:).^2);
    model2_sum = model2_sum + res2.^2;
    SSE = sum(sum(model2_sum));
    % disp(SSE);

    % Save the best model
    if SSE < best_SSE2
        best_SSE2 = SSE;
        best_alpha2 = alpha2;
        best_beta2 = beta2;
        best_y2 = y2;
    end

    aic2 = 2*4 + 19.*log(model2_sum./19);
    bic2 = 4*log(19) + 19.*log(model2_sum./19);

end

%% Fit Model 3: f = alpha + beta * t
% Cross validation loop
best_SSE3 = realmax;
best_alpha3 = zeros(size(s,1), size(s,2));
best_beta3 = zeros(size(s,1), size(s,2));
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
    coeff =pinv(X' * X) * X' * Y;
    alpha3 = coeff(1);
    beta3 = coeff(2);

    
    % There seems to be some NaN values in the alpha3 and beta3 variables,
    % so these need to be set to zero
    alpha3(isnan(alpha3)) = 0;
    beta3(isnan(beta3)) = 0;

    % Calculate sum of squares
    model3_sum = zeros(68);
    res3 = validation_set_targets(:,:) - (alpha3 + beta3.*validation_set(:,:));
    model3_sum = model3_sum + res3.^2;
    model3_sum(isinf(model3_sum)) = 0; % Set inifity values to zeros 
    SSE = sum(sum(model3_sum));
    %disp(SSE);

    % Save the best model
    if SSE < best_SSE3
        best_SSE3 = SSE;
        best_alpha3 = alpha3;
        best_beta3 = beta3;
    end
    
    aic3 = 2*3 + 19.* log(model3_sum./19);
    bic3 = 3*log(19) + 19.*log(model3_sum./19);
end

%% Fit Model 4: f = alpha + beta * t + y * t^2
% Cross validation loop
best_SSE4 = realmax;
best_alpha4 = zeros(size(s,1), size(s,2));
best_beta4 = zeros(size(s,1), size(s,2));
best_y4 = zeros(size(s,1), size(s,2));
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

    % Solve for the alpha and beta coefficients
    coeff = pinv(X' * X) * X' * Y;
    alpha4 = coeff(1);
    beta4 = coeff(2);
    y4 = coeff(3);


    % There seems to be some NaN values in the alpha3 and beta3 variables,
    % so these need to be set to zero
    alpha4(isnan(alpha4)) = 0;
    beta4(isnan(beta4)) = 0;
    y4(isnan(y4)) = 0;
    
    % Calculate sum of square errors
    model4_sum = zeros(68);
    res4 = validation_set_targets(:,:) - (alpha4 + beta4.*validation_set(:,:) + y4.*validation_set(:,:).^2);
    model4_sum = model4_sum + res4.^2;
    model4_sum(isnan(model4_sum)) = 0; % Set NaN values to zeros again 
    model4_sum(isinf(model4_sum)) = 0; % Set inifity values to zeros 
    SSE = sum(sum(model4_sum));
    %disp(SSE);

    % Save the best model
    if SSE < best_SSE4
        best_SSE4 = SSE;
        best_alpha4 = alpha4;
        best_beta4 = beta4;
        best_y4 = y4;
    end

    % AIC/BIC calculation
    aic4 = 2*4 + 19.* log(model4_sum./19);
    bic4 = 4*log(19) + 19.*log(model4_sum./19);
   
end


%% Fit Model 5: f = alpha + beta * s + y * t
% Cross validation loop
best_SSE5 = realmax;
best_alpha5 = zeros(size(s,1), size(s,2));
best_beta5 = zeros(size(s,1), size(s,2));
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
    Y = reshape(training_set_targets(:,:,:), 68 * 68 * 18, 1);

    % Solve for the alpha and beta coefficients
    coeff = pinv(X' * X) * X' * Y;
    alpha5 = coeff(1);
    beta5 = coeff(2);
    y5 = coeff(3);
        
    % There seems to be some NaN values in the alpha3 and beta3 variables,
    % so these need to be set to zero
    alpha5(isnan(alpha5)) = 0;
    beta5(isnan(beta5)) = 0;
    y5(isnan(y5)) = 0;
    
    % Calculate sum of square errors
    model5_sum = zeros(68);
    res5 = validation_set_targets(:,:) - (alpha5 + beta5.*validation_set_s(:,:) + y5.*validation_set_t(:,:));
    model5_sum = model5_sum + res5.^2;
    model5_sum(isnan(model5_sum)) = 0; % Set NaN values to zeros again 
    model5_sum(isinf(model5_sum)) = 0; % Set inifity values to zeros 
    SSE = sum(sum(model5_sum));
    %disp(SSE);
    
    % Save the best model
    if SSE < best_SSE5
        best_SSE5 = SSE;
        best_alpha5 = alpha5;
        best_beta5 = beta5;
    end

    % AIC/BIC
    aic5 = 2*4 + 19.* log(model5_sum./19);
    bic5 = 4*log(19) + 19.*log(model5_sum./19);
end

disp("Model 1 :" + best_SSE1)
disp("Model 2 :" + best_SSE2)
disp("Model 3 :" + best_SSE3)
disp("Model 4 :" + best_SSE4)
disp("Model 5 :" + best_SSE5)










