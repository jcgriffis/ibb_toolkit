function [cv_results] = run_svr_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag)

% Use K-fold cross-validation to tune SVR hyper-parameters
% Joseph Griffis 2023

if cfg.optimize_hyperparams == 1

    % Get hyper-parameter options and cross-validation partitions
    cv_results.hp_opt = get_hp_opts(cfg, length(y_train));
    
    % Fit model with optimization
    alpha = []; % sometimes model returns empty results despite convergence
    while isempty(alpha)
        mdl = fitrsvm(x_train,y_train,'KernelFunction', cfg.kernel,...
            'OptimizeHyperparameters', {'BoxConstraint', 'KernelScale'},...
            'HyperparameterOptimizationOptions', cv_results.hp_opt);
        if isnan(mdl.Bias)
            alpha = [];
        else
            alpha = mdl.Alpha;
        end
    end
    
    % Get parameters
    opt_c = mdl.BoxConstraints(1);
    opt_g = mdl.KernelParameters.Scale;

else
    
    % Fit model with no optimization
    alpha = []; % sometimes model returns empty results despite convergence
    while isempty(alpha)
        mdl = fitrsvm(x_train,y_train, 'KernelFunction', cfg.kernel);
        if isnan(mdl.Bias)
            alpha = [];
        else
            alpha = mdl.Alpha;
        end
    end
    
    % Get parameters
    opt_c = mdl.BoxConstraints(1);
    opt_g = mdl.KernelParameters.Scale;

end

% Get predictions
y_pred = predict(mdl, x_test);
if isfield(cfg, 'Cy')
    y_pred = (y_pred .* cfg.Sy) + cfg.Cy;
    y_test = (y_test .* cfg.Sy) + cfg.Cy;
end

% Store predictions and observations
cv_results.pred_y(cfg.test_set==1,i,j) = y_pred; % Predicted outcome for test set
cv_results.obs_y(cfg.test_set==1,i,j) = y_test; % Observed outcome for test set
cv_results.mse(i,j) = get_mse(y_test, y_pred); % MSE

if perm_flag == 0
    if strcmp(cfg.model_spec, 'kernsvr')
        cv_results.coeff(:,i,j) = (mdl.Alpha.'*mdl.SupportVectors); % Compute betas using sensitivity mapping method (Zhang et al., 2014 - Human Brain Mapping)
    else
        cv_results.coeff(:,i,j) = mdl.Beta;
    end
    cv_results.explained(i,j) = get_explained_variance(y_test, y_pred); % explained variance
    cv_results.r2_ss(i,j) = get_model_r2(y_test, y_pred); % sum-of-squares R-squared
    cv_results.corr(i,j) = corr(y_test, y_pred); % Correlation of predicted and observed    
    cv_results.opt_c(i,j) = opt_c;
    cv_results.opt_g(i,j) = opt_g;
end

end