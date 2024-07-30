function [XS, betas, pred_y, classrate, vip_score, roc_auc, pred_score, opt_thresh] = fitplsda(X, Y, opt_k, cost)

% Fits a PLSR rgression model with specified dimensionality and computes R-squared
% Joseph Griffis 2023
    
% Fit PLS regression
[XL,yl,XS,~,betas,~,~,stats] = plsregress(X, Y, opt_k);

% Get vip scores
W0 = stats.W ./ sqrt(sum(stats.W.^2,1));
p = size(XL,1);
sumSq = sum(XS.^2,1).*sum(yl.^2,1);
vip_score = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));

% Get predicted scores for test set
pred_y = [ones(length(Y), 1) X]*betas; % get fitted Y

% Get ROC AUC 
roc_obj = rocmetrics(Y, pred_y, 1);    
roc_auc = roc_obj.AUC;

% Get optimal threshold
opt_thresh = get_optimal_threshold(Y, roc_obj, cost);

% Store raw predictions
pred_score = pred_y;

% Convert to labels
pred_y(:) = 0;
pred_y(pred_score >= opt_thresh) = 1;
pred_y(pred_score < opt_thresh) = -1;

% Compute accuracy for each group and for full sample
classrate(2) = numel(intersect(find(pred_y==1), find(Y==1)))./numel(find(Y==1));
classrate(1) = numel(intersect(find(pred_y==-1), find(Y==-1)))./numel(find(Y==-1));
classrate(3) = numel(find(pred_y==Y))./numel(find(Y));

end