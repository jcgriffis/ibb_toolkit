function [out] = run_boot_classmdl(X, Y, lambda, learner, reg_type, cost, standardize)

% Run bootstrap analyses
% Joseph Griffis 2024

% Get model params
if learner == 1
    learner = 'svm';
else
    learner = 'logistic';
end
if reg_type == 1
    reg_type = 'lasso';
else
    reg_type = 'ridge';
end

if standardize == 1 || standardize == 3
    [X] = normalize(X);
    X(isnan(X))=0; % Since normalize will cause columns to become NaN if they are all 0
end

% Get costs
S.ClassNames = [1, -1];
S.ClassificationCosts = cost;

% Fit final model
mdl = fitclinear(X,Y,...
    'ObservationsIn', 'rows',...
    'Regularization', reg_type,...
    'Lambda', lambda,...
    'Learner', learner,...
    'BetaTolerance',0,...
    'Cost', S);

% Get betas
betas = mdl.Beta;

% Get final model classification rates
[pred_y, pred_score] = predict(mdl, X, 'ObservationsIn', 'rows');

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
gof(4,1) = roc_obj.AUC;    

out = [gof; betas];

end
