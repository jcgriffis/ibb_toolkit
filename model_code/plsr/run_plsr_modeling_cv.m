function cv_results = run_plsr_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag)

% Run hyperparamter optimization and then fit model 
% Joseph Griffis 2024

%%%% run CV hyperparameter optimmization for number of PLS components
if cfg.optimize_hyperparams == 1
    switch cfg.hp_opt.cv_type 
        case 'LOOCV' % Leave-one-out
            disp('Running LOOCV optimization for PLS components...');
            [opt_k, ~] = run_loocv_plsr(x_train, y_train, cfg); % loo cv to identify opt_k for fixed effects and get out-of-sample prediction
        case 'KFold' % K-Fold 
            disp('Running KFold optimization for PLS component');
            disp(['Number of folds: ' num2str(cfg.hp_opt.folds)]);
            disp(['Number of repeated ' num2str(cfg.hp_opt.folds) '-fold splits: ' num2str(cfg.hp_opt.repeats)]);
            [opt_k, ~] = run_kfold_plsr(x_train, y_train, cfg); % k-fold cv to identify opt_k for fixed effects and get out-of-sample prediction
    end
else
    opt_k = cfg.opt_k;
end
%%%%% Get fitted model with optimal hyperparameters
[~, betas, ~, ~, vip_score] = fitplsrm(x_train, y_train, opt_k);

% Get predicted scores for test set
y_pred = [ones(length(y_test), 1) x_test]*betas; % get fitted Y

% Put predictions and observations back in original units if needed
if isfield(cfg, 'Cy')
    y_pred = (y_pred .* cfg.Sy) + cfg.Cy;
    y_test = (y_test .* cfg.Sy) + cfg.Cy;
end

% Save relevant results in cv_results structure
cv_results.pred_y(cfg.test_set==1,i,j) = y_pred; % Predicted outcome for test set
cv_results.obs_y(cfg.test_set==1,i,j) = y_test; % Observed outcome for test set
cv_results.mse(i,j) = get_mse(y_test, y_pred); % MSE

if perm_flag == 0
    cv_results.explained(i,j) = get_explained_variance(y_test, y_pred); % explained variance
    cv_results.mse(i,j) = get_mse(y_test, y_pred); % MSE
    cv_results.r2_ss(i,j) = get_model_r2(y_test, y_pred); % sum-of-squares R-squared
    cv_results.corr(i,j) = corr(y_test, y_pred); % Correlation of predicted and observed
    cv_results.coeff(:,i,j) = betas(2:end); % betas from training set
    cv_results.beta_0(i,j) = betas(1); % Intercept from training set
    cv_results.vip(:,i,j) = vip_score; % VIP scores from training set
    cv_results.optk(i,j) = opt_k; % optimal component number from training set
end
disp(['Outer Fold Test Set (' num2str(i) ',' num2str(j) ')'])                

end

