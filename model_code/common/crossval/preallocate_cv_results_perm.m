function cv_results = preallocate_cv_results_perm(row_x, cfg)

% Preallocate relevant fields to store CV results
% Joseph Griffis 2024

disp('Preallocating CV results')

% CV Outputs
if strcmp(cfg.cv.cv_type, 'LOOCV')
    cfg.cv.folds = row_x;
    cfg.cv.repeats = 1;
end
   
% Confound regression outputs
if isfield(cfg, 'confounds') && ~isempty(cfg.confounds) % Confound regression
    cv_results.r2_confound = zeros(cfg.cv.repeats, cfg.cv.folds);
    cv_results.coeff_confound = zeros(size(cfg.confounds,2)+1,cfg.cv.repeats, cfg.cv.folds);
end

% Common outputs
cv_results.all.pred_y = zeros(row_x, cfg.cv.repeats).*NaN;
cv_results.all.obs_y = zeros(row_x, cfg.cv.repeats).*NaN;
cv_results.pred_y = zeros(row_x, cfg.cv.repeats, cfg.cv.folds).*NaN;
cv_results.obs_y = zeros(row_x, cfg.cv.repeats, cfg.cv.folds).*NaN;

% Outputs for continuous outcomes
if cfg.cat_Y == 0
    
    % Results for full sample (collapsed across folds)
    cv_results.mse = zeros(cfg.cv.repeats, cfg.cv.folds);
    cv_results.avg.mse = zeros(cfg.cv.repeats, 1);
    cv_results.pred_score = zeros(row_x, cfg.cv.repeats, cfg.cv.folds).*NaN;

else

    % Results for full sample (collapsed across folds)    
    cv_results.classrate1 = zeros(cfg.cv.repeats, cfg.cv.folds);
    cv_results.classrate2 = zeros(cfg.cv.repeats, cfg.cv.folds); 
    cv_results.roc_auc = zeros(cfg.cv.repeats, cfg.cv.folds);
    cv_results.avg.classrate1 = zeros(cfg.cv.repeats, 3);
    cv_results.avg.classrate2 = zeros(cfg.cv.repeats, 3);    
    cv_results.avg.roc_auc = zeros(cfg.cv.repeats, 1);

end