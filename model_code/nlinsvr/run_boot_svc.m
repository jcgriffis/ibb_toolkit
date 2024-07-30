function [out] = run_boot_svc(X, Y, C, gamma, kernel, standardize, cost)

% Run bootstrap analyses
% Joseph Griffis 2024

% Set kernel function
if kernel == 0
    kfun = 'linear';
else
    kfun = 'rbf';
end

if standardize == 1 || standardize == 3
    [X] = normalize(X);
    X(isnan(X))=0; % Since normalize will cause columns to become NaN if they are all 0
end

S.ClassNames = [1, -1];
S.ClassificationCosts = cost;

alpha = [];
while isempty(alpha)
    % Fit final model with cross-validation optimized hyper-parameters
    mdl = fitcsvm(X,Y,'KernelFunction', kfun,...
    'BoxConstraint', C,...
    'KernelScale', gamma,...
    'Cost', S);
    if isnan(mdl.Bias)
        alpha = [];
    else
        alpha = mdl.Alpha;
    end
end

% Compute betas using sensitivity mapping method (Zhang et al., 2014 - Human Brain Mapping)
if strcmp(kfun, 'rbf')
    betas(:,1) = (mdl.Alpha.'*(mdl.SupportVectors.*mdl.SupportVectorLabels));
else
    betas(:,1) = mdl.Beta;
end

% Get final model classification rates
[pred_y, pred_score] = predict(mdl, X);

% Get ROC AUC 
roc_obj = rocmetrics(Y, pred_score(:,2), 1);    

% Get optimal threshold
opt_thresh = get_optimal_threshold(Y, roc_obj, cost);

% Assign labels using optimized threshold
pred_y(:) = 0;
pred_y(pred_score(:,2) >= opt_thresh) = 1;
pred_y(pred_score(:,2) < opt_thresh) = -1;

% Compute accuracy for each group and for full sample
gof(2,1) = numel(intersect(find(pred_y==1), find(Y==1)))./numel(find(Y==1));
gof(1,1) = numel(intersect(find(pred_y==-1), find(Y==-1)))./numel(find(Y==-1));
gof(3,1) = numel(find(pred_y==Y))./numel(find(Y));
   
% Get ROC AUC 
roc_obj = rocmetrics(Y, pred_score(:,2), 1);    
gof(4,1) = roc_obj.AUC;    
    
out = [gof; betas];

end