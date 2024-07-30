function model_results = run_plsr_modeling(X, Y, cfg, model_results)

% Run hyperparamter optimization and then fit model 
% Joseph Griffis 2024
disp(cfg.hp_opt.cv_type);

%%%% run CV hyperparameter optimmization for number of PLS components
if cfg.optimize_hyperparams == 1
    switch cfg.hp_opt.cv_type
        case 'LOOCV' % Leave-one-out
            disp('Running LOOCV optimization for PLS components...');
            [model_results.opt_k, model_results.hp_opt_mse] = run_loocv_plsr(X, Y, cfg); % loo cv to identify opt_k for fixed effects and get out-of-sample prediction
            disp(['The optimal number of PLS components is ' num2str(model_results.opt_k) '. Fitting full model with optimal components.']);
        case 'KFold' % K-Fold 
            disp('Running KFold optimization for PLS component');
            disp(['Number of folds: ' num2str(cfg.hp_opt.folds)]);
            disp(['Number of repeated ' num2str(cfg.hp_opt.folds) '-fold splits: ' num2str(cfg.hp_opt.repeats)]);
            [model_results.opt_k, model_results.hp_opt_mse] = run_kfold_plsr(X, Y, cfg); % k-fold cv to identify opt_k for fixed effects and get out-of-sample prediction
            disp(['The optimal number of PLS components is ' num2str(model_results.opt_k) '. Fitting full model with optimal components.']);
    end
else
    model_results.opt_k = cfg.opt_k;
end
   
%%%%% Get fitted model with optimal hyperparameters
[model_results.XS, model_results.coeff, model_results.pred_y, ~, model_results.vip_score, ~] = fitplsrm(X, Y, model_results.opt_k);

% Put predictions and observations back in original units if needed
if isfield(model_results, 'Sy')
    model_results.pred_y = (model_results.pred_y .* model_results.Sy) + model_results.Cy;
    model_results.obs_y = (Y .* model_results.Sy) + model_results.Cy;
else
    model_results.obs_y = Y;
end

% Get performance metrics
model_results.r2 = get_model_r2(model_results.obs_y, model_results.pred_y);
model_results.mse = get_mse(model_results.obs_y, model_results.pred_y);
disp(['The full model explains ' num2str(100*model_results.r2) ' % of the variance in the response.']);

% Get betas
model_results.beta_0 = model_results.coeff(1);
model_results.coeff = model_results.coeff(2:end);

end

