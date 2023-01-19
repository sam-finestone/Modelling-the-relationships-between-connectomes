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
% Calucluate indirect structural connectivity matrix
t = zeros(size(s));
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

%% Additional model 9 f = alpha + beta * t + gamma * s.^2
alpha9 = zeros(size(s,1),size(s,2)); 
beta9 = zeros(size(s,1),size(s,2)); 
gamma9 = zeros(size(s,1),size(s,2)); 

for i = 1:size(s, 1)
    for j = 1:size(s, 2)
        s_slice = reshape(s(i,j,:), 19, 1);
        t_slice = reshape(t(i,j,:), 19, 1);
        X = [ones(19, 1) t_slice (s_slice.^2)];
        Y = reshape(f(i,j,:), 19, 1);

        coeff = inv(X' * X) * X' * Y;
        alpha9(i,j) = coeff(1);
        beta9(i,j) = coeff(2);
        gamma9(i,j) = coeff(3);
    end
end

SSE9 = zeros(19,1);
for k = 1:19
    model9_sum = zeros(68);
    res9 = f(:,:,k) - (alpha9 + beta9.*t(:,:,k) + gamma9.*s(:,:,k).^2);
    model9_sum = model9_sum + res9.^2;
    model9_sum(isinf(model9_sum)) = 0;
    model9_sum(isnan(model9_sum)) = 0;
    SSE9(k) = sum(sum(model9_sum));
end

best_SSE9 = realmax;
best_alpha9 = zeros(size(s,1), size(s,2));
best_beta9 = zeros(size(s,1), size(s,2));
best_gamma9 = zeros(size(s,1), size(s,2));
best_aic9 = 0;
best_bic9 = 0;

for cv = 1:19
    f_test = f(:,:,cv);
    f_train = f(:,:,[1:cv-1,cv+1:19]);
    s_test = s(:,:,cv);
    s_train = s(:,:,[1:cv-1,cv+1:19]);
    t_test = t(:,:,cv);
    t_train = t(:,:,[1:cv-1,cv+1:19]);
    
    cv_alpha9 = zeros(size(s,1),size(s,2));
    cv_beta9 = zeros(size(s,1),size(s,2));
    cv_gamma9 = zeros(size(s,1),size(s,2));
    
    for i = 1:size(s, 1)
        for j = 1:size(s, 2)
            s_slice = reshape(s_train(i,j,:), 18, 1);
            t_slice = reshape(t_train(i,j,:), 18, 1);
            cv_X = [ones(18, 1) t_slice (s_slice.^2)];
            cv_Y = reshape(f_train(i,j,:), 18, 1);

            cv_coeff = inv(cv_X' * cv_X) * cv_X' * cv_Y;
            cv_alpha9(i,j) = cv_coeff(1);
            cv_beta9(i,j) = cv_coeff(2); 
            cv_gamma9(i,j) = cv_coeff(3);

        end
    end
    
    cv_model9_sum = zeros(68);
    cv_res9 = f_test - (cv_alpha9 + cv_beta9.*t_test + cv_gamma9.*s_test.^2);
    cv_model9_sum = cv_model9_sum + cv_res9.^2;
    cv_model9_sum(isinf(cv_model9_sum)) = 0;
    cv_model9_sum(isnan(cv_model9_sum)) = 0;
    cv_SSE = sum(cv_model9_sum,'all');
    
    aic9 = 2*4 + 19.* log(cv_model9_sum./19);
    bic9 = 4*log(19) + 19.*log(cv_model9_sum./19);
    
    if cv_SSE < best_SSE9
        best_SSE9 = cv_SSE;
        best_alpha9 = cv_alpha9;
        best_beta9 = cv_beta9;
        best_gamma9 = cv_gamma9;
        best_aic9 = MeanExcludingInf(aic9);
        best_bic9 = MeanExcludingInf(bic9);
    end
end  

%% Additional model 6 f = alpha + beta * t.^3
alpha6 = zeros(size(s,1),size(s,2)); 
beta6 = zeros(size(s,1),size(s,2)); 

for i = 1:size(s, 1)
    for j = 1:size(s, 2)
        t_slice = reshape(s(i,j,:), 19, 1);
        X = [ones(19, 1) (t_slice.^3)];
        Y = reshape(f(i,j,:), 19, 1);

        coeff = pinv(X' * X) * X' * Y;
        alpha6(i,j) = coeff(1);
        beta6(i,j) = coeff(2);
    end
end

SSE6 = zeros(19,1);
for k = 1:19
    model6_sum = zeros(68);
    res6 = f(:,:,k) - (alpha6 + beta6.*t(:,:,k).^3);
    model6_sum = model6_sum + res6.^2;
    model6_sum(isinf(model6_sum)) = 0;
    model6_sum(isnan(model6_sum)) = 0;
    SSE6(k) = sum(sum(model6_sum));
end

aic6 = 2*3 + 19.*log(model6_sum./19);
bic6 = 3*log(19) + 19.*log(model6_sum./19);

cv_SSE6 = zeros(19,1);
for cv = 1:19
    f_test = f(:,:,cv);
    f_train = f(:,:,[1:cv-1,cv+1:19]);
    t_test = t(:,:,cv);
    t_train = t(:,:,[1:cv-1,cv+1:19]);
    
    cv_alpha6 = zeros(size(s,1),size(s,2));
    cv_beta6 = zeros(size(s,1),size(s,2));
    
    for i = 1:size(s, 1)
        for j = 1:size(s, 2)
            t_slice = reshape(t_train(i,j,:), 18, 1);
            cv_X = [ones(18, 1) (t_slice.^3)];
            cv_Y = reshape(f_train(i,j,:), 18, 1);

            cv_coeff = inv(cv_X' * cv_X) * cv_X' * cv_Y;
            cv_alpha6(i,j) = cv_coeff(1);
            cv_beta6(i,j) = cv_coeff(2); 

        end
    end
    
    cv_model6_sum = zeros(68);
    cv_res6 = f_test - (cv_alpha6 + cv_beta6.*t_test.^3);
    cv_model6_sum = cv_model6_sum + cv_res6.^2;
    cv_model6_sum(isinf(cv_model6_sum)) = 0;
    cv_model6_sum(isnan(cv_model6_sum)) = 0;
    cv_SSE6(cv) = sum(cv_model6_sum,'all');
end    

% plot the relationship
figure; 
for i=1:10   
    position = randi([1,68],2,1);
    fval = [];
    tval = [];
    for k = 1:19
        fval = [fval,f(position(1),position(2),k)];
        tval = [tval,t(position(1),position(2),k)];
    end
    
    subplot(2,5,i);
    x = linspace(0,0.7);
    a = alpha6(position(1),position(2));
    b = beta6(position(1),position(2));
    plot(x,a+b.*x.^3);
    hold on;
    
    scatter(tval,fval);
    title(sprintf('%.0f,%.0f', position(1),position(2)));
end

%% Additional model 7 f = alpha + beta * s.^2

alpha7 = zeros(size(s,1),size(s,2)); 
beta7 = zeros(size(s,1),size(s,2)); 

for i = 1:size(s, 1)
    for j = 1:size(s, 2)
        s_slice = reshape(s(i,j,:), 19, 1);
        X = [ones(19, 1) (s_slice.^2)];
        Y = reshape(f(i,j,:), 19, 1);

        coeff = pinv(X' * X) * X' * Y;
        alpha7(i,j) = coeff(1);
        beta7(i,j) = coeff(2);
    end
end

SSE7 = zeros(19,1);
for k = 1:19
    model7_sum = zeros(68);
    res7 = f(:,:,k) - (alpha7 + beta7.*s(:,:,k).^2);
    model7_sum = model7_sum + res7.^2;
    model7_sum(isinf(model7_sum)) = 0;
    model7_sum(isnan(model7_sum)) = 0;
    SSE7(k) = sum(sum(model7_sum));
end
aic7 = 2*3 + 19.*log(model7_sum./19);
bic7 = 3*log(19) + 19.*log(model7_sum./19);


cv_SSE7 = zeros(19,1);
for cv = 1:19
    f_test = f(:,:,cv);
    f_train = f(:,:,[1:cv-1,cv+1:19]);
    s_test = s(:,:,cv);
    s_train = s(:,:,[1:cv-1,cv+1:19]);
    
    cv_alpha7 = zeros(size(s,1),size(s,2));
    cv_beta7 = zeros(size(s,1),size(s,2));
    
    for i = 1:size(s, 1)
        for j = 1:size(s, 2)
            s_slice = reshape(s_train(i,j,:), 18, 1);
            cv_X = [ones(18, 1) (s_slice.^2)];
            cv_Y = reshape(f_train(i,j,:), 18, 1);

            cv_coeff = pinv(cv_X' * cv_X) * cv_X' * cv_Y;
            cv_alpha7(i,j) = cv_coeff(1);
            cv_beta7(i,j) = cv_coeff(2); 

        end
    end
    
    cv_model7_sum = zeros(68);
    cv_res7 = f_test - (cv_alpha7 + cv_beta7.*s_test.^2);
    cv_model7_sum = cv_model7_sum + cv_res7.^2;
    cv_model7_sum(isinf(cv_model7_sum)) = 0;
    cv_model7_sum(isnan(cv_model7_sum)) = 0;
    cv_SSE7(cv) = sum(cv_model7_sum,'all');
end    


% plot the relationship
figure; 
for i=1:10   
    position = randi([1,68],2,1);
    fval = [];
    sval = [];
    for k = 1:19
        fval = [fval,f(position(1),position(2),k)];
        sval = [sval,s(position(1),position(2),k)];
    end
    
    subplot(2,5,i);
    x = linspace(0,0.7);
    a = alpha7(position(1),position(2));
    b = beta7(position(1),position(2));
    plot(x,a+b.*x.^2);
    hold on;
    
    scatter(sval,fval);
    title(sprintf('(%.0f , %.0f)', position(1),position(2)));
end

%% Additional model 8 f = alpha + beta * s + gamma * t.^2

alpha8 = zeros(size(s,1),size(s,2)); 
beta8 = zeros(size(s,1),size(s,2)); 
gamma8 = zeros(size(s,1),size(s,2)); 

for i = 1:size(s, 1)
    for j = 1:size(s, 2)
        s_slice = reshape(s(i,j,:), 19, 1);
        t_slice = reshape(t(i,j,:), 19, 1);
        X = [ones(19, 1) s_slice (t_slice.^2)];
        Y = reshape(f(i,j,:), 19, 1);

        coeff = inv(X' * X) * X' * Y;
        alpha8(i,j) = coeff(1);
        beta8(i,j) = coeff(2);
        gamma8(i,j) = coeff(3);
    end
end

SSE8 = zeros(19,1);
for k = 1:19
    model8_sum = zeros(68);
    res8 = f(:,:,k) - (alpha8 + beta8.*s(:,:,k) + gamma8.*t(:,:,k).^2);
    model8_sum = model8_sum + res8.^2;
    model8_sum(isinf(model8_sum)) = 0;
    model8_sum(isnan(model8_sum)) = 0;
    SSE8(k) = sum(sum(model8_sum));
end

cv_SSE8 = zeros(19,1);
for cv = 1:19
    f_test = f(:,:,cv);
    f_train = f(:,:,[1:cv-1,cv+1:19]);
    s_test = s(:,:,cv);
    s_train = s(:,:,[1:cv-1,cv+1:19]);
    t_test = t(:,:,cv);
    t_train = t(:,:,[1:cv-1,cv+1:19]);
    
    cv_alpha8 = zeros(size(s,1),size(s,2));
    cv_beta8 = zeros(size(s,1),size(s,2));
    cv_gamma8 = zeros(size(s,1),size(s,2));
    
    for i = 1:size(s, 1)
        for j = 1:size(s, 2)
            s_slice = reshape(s_train(i,j,:), 18, 1);
            t_slice = reshape(t_train(i,j,:), 18, 1);
            cv_X = [ones(18, 1) s_slice (t_slice.^2)];
            cv_Y = reshape(f_train(i,j,:), 18, 1);

            cv_coeff = inv(cv_X' * cv_X) * cv_X' * cv_Y;
            cv_alpha8(i,j) = cv_coeff(1);
            cv_beta8(i,j) = cv_coeff(2); 
            cv_gamma8(i,j) = cv_coeff(3);

        end
    end
    
    cv_model8_sum = zeros(68);
    cv_res8 = f_test - (cv_alpha8 + cv_beta8.*s_test + cv_gamma8.*t_test.^2);
    cv_model8_sum = cv_model8_sum + cv_res8.^2;
    cv_model8_sum(isinf(cv_model8_sum)) = 0;
    cv_model8_sum(isnan(cv_model8_sum)) = 0;
    cv_SSE8(cv) = sum(cv_model8_sum,'all');
end    

%%
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