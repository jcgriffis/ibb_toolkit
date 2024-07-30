function [cv_results] = run_svc_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag)

% Use K-fold cross-validation to tune SVC hyper-parameters and fit model in nested CV scheme
% Joseph Griffis 2023

% Get costs
S.ClassNames = [-1, 1];
S.ClassificationCosts = cfg.cost;

if cfg.optimize_hyperparams == 1

    % Get hyper-parameter options and cross-validation partitions
    cv_results.hp_opt = get_hp_opts(cfg, length(y_train));
    
    % Fit model with optimization
    alpha = []; % sometimes model returns empty results despite convergence
    while isempty(alpha)
        mdl = fitcsvm(x_train,y_train,'KernelFunction', cfg.kernel,...
            'OptimizeHyperparameters', {'BoxConstraint', 'KernelScale'},...
            'HyperparameterOptimizationOptions', cv_results.hp_opt, ...
            'Cost', S);
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
        mdl = fitcsvm(x_train,y_train, 'KernelFunction', cfg.kernel,...
            'Cost', S);
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

% Get optimal threshold from training set
[~, pred_score] = predict(mdl, x_train);
roc_obj = rocmetrics(y_train, pred_score(:,2), 1);    
opt_thresh = get_optimal_threshold(y_train, roc_obj, cfg.cost);

% Get predictions and coefficients
[pred_y, pred_score] = predict(mdl, x_test);

% Assign labels using optimized threshold
pred_y(:) = 0;
pred_y(pred_score(:,2) >= opt_thresh) = 1;
pred_y(pred_score(:,2) < opt_thresh) = -1;

% Store predicted and observed
cv_results.pred_y(cfg.test_set==1,i,j) = pred_y; % Predicted outcome for test set
cv_results.obs_y(cfg.test_set==1,i,j) = y_test; % Observed outcome for test set
roc_obj = rocmetrics(y_test, pred_score(:,2), 1);    
cv_results.roc_auc(i,j) = roc_obj.AUC;   

if perm_flag == 0
    if strcmp(cfg.model_spec, 'kernsvc')
        cv_results.coeff(:,i,j) = (mdl.Alpha.'*(mdl.SupportVectors.*mdl.SupportVectorLabels)); % Compute betas using sensitivity mapping method (Zhang et al., 2014 - Human Brain Mapping)
    else
        cv_results.coeff(:,i,j) = mdl.Beta;
    end    
    cv_results.pred_score(cfg.test_set==1,i,j) = pred_score(:,2);   
    cv_results.classrate(i,j,2) = numel(intersect(find(pred_y==1), find(y_test==1)))./numel(find(y_test==1));
    cv_results.classrate(i,j,1) = numel(intersect(find(pred_y==-1), find(y_test==-1)))./numel(find(y_test==-1));    
    cv_results.classrate(i,j,3) = numel(find(pred_y==y_test))./numel(find(y_test));
    cv_results.opt_c(i,j) = opt_c;
    cv_results.opt_g(i,j) = opt_g;
end

end