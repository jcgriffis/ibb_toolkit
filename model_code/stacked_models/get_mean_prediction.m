function stacked_results = get_mean_prediction(x_test,y_test, cfg, stacked_results,i,j)

% This function simply averages the predictions across models and returns the result as the predicted outcome
% Joseph Griffis, 2024

% Get average prediction
y_pred = nanmean(x_test, 2);

% Get model performance measures for regression
stacked_results.explained(i,j) = get_explained_variance(y_test, y_pred);
stacked_results.mse(i,j) = get_mse(y_test, y_pred);
stacked_results.corr(i,j) = corr(y_test, y_pred);
stacked_results.r2_ss(i,j) = get_model_r2(y_test, y_pred);
stacked_results.pred_y(cfg.test_set==1,i,j) = y_pred;
stacked_results.obs_y(cfg.test_set==1,i,j) = y_test;

end