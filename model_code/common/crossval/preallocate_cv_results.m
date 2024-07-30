function cv_results = preallocate_cv_results(row_x, col_x, cfg)

% Preallocate relevant fields to store CV results
% Joseph Griffis 2024

disp('Preallocating CV results')

% CV Outputs
if strcmp(cfg.cv.cv_type, 'LOOCV')
    cfg.cv.folds = row_x;
    cfg.cv.repeats = 1;
end
cv_results.train_set = zeros(row_x, cfg.cv.repeats, cfg.cv.folds);
cv_results.test_set = zeros(row_x, cfg.cv.repeats, cfg.cv.folds);

% Model specific outputs
switch cfg.model_spec
    
    case {'plsr', 'pls_da'}
    
        % PLS-specific model outputs
        cv_results.beta_0 = zeros(cfg.cv.repeats, cfg.cv.folds);
        cv_results.vip = zeros(col_x, cfg.cv.repeats, cfg.cv.folds);
        cv_results.optk = zeros(cfg.cv.repeats, cfg.cv.folds);

    case {'ridge', 'lasso', 'rlinsvr', 'logistic_ridge', 'logistic_lasso', 'rlinsvc'} 
        
        % RegLM-specific model outputs
        cv_results.beta_0 = zeros(cfg.cv.repeats, cfg.cv.folds);
        cv_results.lambda = zeros(cfg.cv.repeats, cfg.cv.folds);

    case {'linsvr', 'kernsvr'}
        
        % SVR-specific outputs
        cv_results.C = zeros(cfg.cv.repeats, cfg.cv.folds);
        cv_results.epsilon = zeros(cfg.cv.repeats, cfg.cv.folds);
        cv_results.gamma = zeros(cfg.cv.repeats, cfg.cv.folds);
        cv_results.alpha = zeros(cfg.cv.repeats, cfg.cv.folds);

    case {'linsvc', 'kernsvc'}

        % SVC-specific outputs
        cv_results.C = zeros(cfg.cv.repeats, cfg.cv.folds);
        cv_results.gamma = zeros(cfg.cv.repeats, cfg.cv.folds);
        cv_results.alpha = zeros(cfg.cv.repeats, cfg.cv.folds);  

    case {'censemble', 'rensemble'}

        % Ensemble specific outputs
        cv_results.method = cell(cfg.cv.repeats, cfg.cv.folds);
        cv_results.n_learn = zeros(cfg.cv.repeats, cfg.cv.folds);

end
   
% Confound regression outputs
if isfield(cfg, 'confounds') && ~isempty(cfg.confounds) % Confound regression
    cv_results.avg.r2_confound = zeros(cfg.cv.repeats, 1);
    cv_results.r2_confound = zeros(cfg.cv.repeats, cfg.cv.folds);
    cv_results.coeff_confound = zeros(size(cfg.confounds,2)+1,cfg.cv.repeats, cfg.cv.folds);
end

% Coefficients for other models are indexed differently due to interactions in OLSR
if ~strcmp(cfg.model_spec, 'olsr')
    cv_results.avg.coeff = zeros(col_x, cfg.cv.repeats).*NaN;
    cv_results.coeff = zeros(col_x, cfg.cv.repeats, cfg.cv.folds);   
else
        if isfield(cfg, 'model_dim')
            cv_results.coeff = zeros(cfg.mdl_dim, cfg.cv.repeats, cfg.cv.folds).*NaN;
            cv_results.avg.coeff = zeros(cfg.mdl_dim, cfg.cv.repeats).*NaN;
        else
            cv_results.coeff = zeros(col_x, cfg.cv.repeats, cfg.cv.folds).*NaN;
            cv_results.avg.coeff = zeros(col_x, cfg.cv.repeats).*NaN;
        end     
end

% Common outputs
cv_results.all.pred_y = zeros(row_x, cfg.cv.repeats).*NaN;
cv_results.all.obs_y = zeros(row_x, cfg.cv.repeats).*NaN;
cv_results.pred_y = zeros(row_x, cfg.cv.repeats, cfg.cv.folds).*NaN;
cv_results.obs_y = zeros(row_x, cfg.cv.repeats, cfg.cv.folds).*NaN;
if cfg.cat_Y == 1
    cv_results.all.pred_score = zeros(row_x, cfg.cv.repeats).*NaN;    
    cv_results.pred_score = zeros(row_x, cfg.cv.repeats, cfg.cv.folds).*NaN;
end

% Outputs for continuous outcomes
if cfg.cat_Y == 0
    
    % Results for full sample (collapsed across folds)
    cv_results.avg.explained = zeros(cfg.cv.repeats, 1);
    cv_results.avg.r2_ss = zeros(cfg.cv.repeats, 1);
    cv_results.avg.mse = zeros(cfg.cv.repeats, 1);
    cv_results.avg.corr = zeros(cfg.cv.repeats, 1);

    % Fold-level results
    cv_results.explained = zeros(cfg.cv.repeats, cfg.cv.folds);
    cv_results.r2_ss = zeros(cfg.cv.repeats, cfg.cv.folds);
    cv_results.mse = zeros(cfg.cv.repeats, cfg.cv.folds);
    cv_results.corr = zeros(cfg.cv.repeats, cfg.cv.folds);

else

    % Results for full sample (collapsed across folds)    
    cv_results.avg.classrate = zeros(cfg.cv.repeats, 3);
    cv_results.avg.roc_auc = zeros(cfg.cv.repeats, 1);

    % Fold-level results  
    cv_results.classrate = zeros(cfg.cv.repeats, cfg.cv.folds, 3);
    cv_results.roc_auc = zeros(cfg.cv.repeats, cfg.cv.folds);   
    cv_results.opt_thresh = zeros(cfg.cv.repeats, cfg.cv.folds);

end