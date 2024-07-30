function cv_results = run_plsda_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag)

% Run hyperparamter optimization and then fit model 
% Joseph Griffis 2024

%%%% run CV hyperparameter optimmization for number of PLS components
if cfg.optimize_hyperparams == 1
    switch cfg.hp_opt.cv_type 
        case 'LOOCV' % Leave-one-out
            disp('Running LOOCV optimization for PLS components...');
            [opt_k, cv_mse] = run_loocv_plsr(x_train, y_train, cfg); % loo cv to identify opt_k for fixed effects and get out-of-sample prediction
        case 'KFold' % K-Fold 
            disp('Running KFold optimization for PLS component');
            disp(['Number of folds: ' num2str(cfg.hp_opt.folds)]);
            disp(['Number of repeated ' num2str(cfg.hp_opt.folds) '-fold splits: ' num2str(cfg.hp_opt.repeats)]);
            [opt_k, cv_mse] = run_kfold_plsr(x_train, y_train, cfg); % k-fold cv to identify opt_k for fixed effects and get out-of-sample prediction
    end
else
    opt_k = 1;
end

%%%%% Get fitted model with optimal hyperparameters
[~, betas, ~, ~, vip_score, ~, ~, opt_thresh] = fitplsda(x_train, y_train, opt_k, cfg.cost);

% Get predicted scores for test set
pred_y = [ones(length(y_test), 1) x_test]*betas; % get fitted Y

%%% Compute accuracy for each group and for full sample

% Get labels using optimal training set threshold
y_label = zeros(length(pred_y),1);
y_label(pred_y >= opt_thresh) = 1; 
y_label(pred_y < opt_thresh) = -1;

% Store predicted and observed labels
cv_results.pred_y(cfg.test_set==1,i,j) = y_label; % Predicted outcome for test set
cv_results.obs_y(cfg.test_set==1,i,j) = y_test; % Observed outcome for test set

% Compute ROC curve on test set
roc_obj = rocmetrics(y_test, pred_y, 1);    
cv_results.roc_auc(i,j) = roc_obj.AUC;   

% Save relevant results in cv_results structure
if perm_flag == 0 
    cv_results.pred_score(cfg.test_set==1,i,j) = pred_y;    
    cv_results.opt_thresh(i,j) = opt_thresh;
    cv_results.classrate(i,j,2) = numel(intersect(find(y_label==1), find(y_test==1)))./numel(find(y_test==1));
    cv_results.classrate(i,j,1) = numel(intersect(find(y_label==-1), find(y_test==-1)))./numel(find(y_test==-1));    
    cv_results.classrate(i,j,3) = numel(find(y_label==y_test))./numel(find(y_test));   
    cv_results.coeff(:,i,j) = betas(2:end); % betas from training set
    cv_results.beta_0(i,j) = betas(1); % intercept from training set
    cv_results.vip(:,i,j) = vip_score; % VIP scores from training set
    cv_results.optk(i,j) = opt_k; % optimal component number from training set
    cv_results.hopt_mse(:,i,j) = cv_mse; % component-wise MSE for training set
end

disp(['Outer Fold Test Set (' num2str(i) ',' num2str(j) ')'])                        

end