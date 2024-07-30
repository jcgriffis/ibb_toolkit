function cv_results = run_ols_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag)
    
% fit OLS regression model and get predictions on test data - limited functionality due to cfg.model_terms not being supported fully
% Joseph Griffis 2024

% Create empty field if it doesn't exist
if ~isfield(cfg, 'categorical_X')
    cfg.categorical_X = [];
end

% Fit model to training dataset
if ~isempty(cfg.categorical_X)
    if isfield(cfg, 'model_terms')
        mdl = fitlm(x_train, y_train, cfg.model_terms, "CategoricalVars", cfg.categorical_X);
    else
        mdl = fitlm(x_train, y_train, "CategoricalVars", cfg.categorical_X);
    end
else
    if isfield(cfg, 'model_terms')
        mdl = fitlm(x_train, y_train, cfg.model_terms);
    else
        mdl = fitlm(x_train, y_train);
    end
end

% Get predicted scores for test set
y_pred = predict(mdl, x_test); % get fitted Y
y_pred(isnan(y_pred)) = mean(y_train); % if there is a NaN, replace with mean of training set

% Put predictions and observations back in original units if needed
if cfg.standardize > 0 
    y_pred = (y_pred .* cfg.Sy) + cfg.Cy;
    y_test = (y_test .* cfg.Sy) + cfg.Cy;
end

% Save relevant results in cv_results structure
cv_results.pred_y(cfg.test_set==1,i,j) = y_pred; % Predicted outcome for test set
cv_results.obs_y(cfg.test_set==1,i,j) = y_test; % Observed outcome for test set
if perm_flag == 0
    cv_results.explained(i,j) = get_explained_variance(y_test, y_pred); % explained variance
    cv_results.mse(i,j) = get_mse(y_test, y_pred); % MSE
    cv_results.r2_ss(i,j) = get_model_r2(y_test, y_pred); % sum-of-squares R-squared
    cv_results.corr(i,j) = corr(y_test, y_pred); % Correlation of predicted and observed            
end

end
