function model_results = run_plsda_modeling(X, Y, cfg, model_results)

% Run hyperparamter optimization and then fit model 
% Joseph Griffis 2024

%%%% run CV hyperparameter optimmization for number of PLS components
if cfg.optimize_hyperparams == 1
    switch cfg.hp_opt.cv_type
        case 'LOOCV' % Leave-one-out
            disp('Running LOOCV optimization for PLS components...');
            [model_results.opt_k, model_results.cv_mse] = run_loocv_plsr(X, Y, cfg); % loo cv to identify opt_k for fixed effects and get out-of-sample prediction
            disp(['The optimal number of PLS components is ' num2str(model_results.opt_k) '. Fitting full model with optimal components.']);
        case 'KFold' % K-Fold 
            disp('Running KFold optimization for PLS component');
            disp(['Number of folds: ' num2str(cfg.hp_opt.folds)]);
            disp(['Number of repeated ' num2str(cfg.hp_opt.folds) '-fold splits: ' num2str(cfg.hp_opt.repeats)]);
            [model_results.opt_k, model_results.cv_mse] = run_kfold_plsr(X, Y, cfg); % k-fold cv to identify opt_k for fixed effects and get out-of-sample prediction
            disp(['The optimal number of PLS components is ' num2str(model_results.opt_k) '. Fitting full model with optimal components.']);
    end
else
    model_results.opt_k = 1;
end
%%%%% Get fitted model with optimal hyperparameters
[model_results.XS, model_results.coeff, model_results.pred_y, model_results.classrate, model_results.vip_score, model_results.roc_auc, model_results.pred_score, model_results.opt_score_thresh] = fitplsda(X, Y, model_results.opt_k, cfg.cost);

% Display main model results
disp(['Classification rate for group -1: ' num2str(round(100.*model_results.classrate(1)),2) '%'])                
disp(['Classification rate for group 1: ' num2str(round(100.*model_results.classrate(2)),2) '% '])                
disp(['The full model has an overall classification rate of ' num2str(round(100.*model_results.classrate(3)),2) '%'])                
disp(['The area under the ROC curve is: ' num2str(model_results.roc_auc)]);

% Add/modify relevant fields in model_results structure
model_results.obs_y = Y;
model_results.beta_0 = model_results.coeff(1);
model_results.coeff = model_results.coeff(2:end);

end
