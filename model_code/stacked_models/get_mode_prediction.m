function stacked_results = get_mode_prediction(x_train, y_train, x_test,y_test, cfg, stacked_results,i,j)

% This function simply averages the predictions across models and returns the result as the predicted outcome
% Joseph Griffis, 2024

% Get optimal threshold from training set
pred_score = mean(x_train,2);
roc_obj = rocmetrics(y_train, pred_score, 1);    
opt_thresh = get_optimal_threshold(y_train, roc_obj, cfg.cost);

% Assign labels using optimized threshold
pred_score = mean(x_test, 2);
pred_y(:) = 0;
pred_y(pred_score >= opt_thresh) = 1;
pred_y(pred_score < opt_thresh) = -1;

% Get model performance measures for regression
stacked_results.classrate(i,j,1) = numel(intersect(find(pred_y==1), find(y_test==1)))./numel(find(y_test==1));
stacked_results.classrate(i,j,2) = numel(intersect(find(pred_y==-1), find(y_test==-1)))./numel(find(y_test==-1));
stacked_results.classrate(i,j,3) = numel(find(pred_y==y_test))./numel(find(y_test));

roc_obj = rocmetrics(y_test, pred_score, 1);    
stacked_results.roc_auc(i,j) = roc_obj.AUC;   
stacked_results.pred_y(cfg.test_set==1,i,j) = pred_y;
stacked_results.pred_score(cfg.test_set==1,i,j) = pred_score; 

end