function cv_results = run_nested_cv(X, Y, cfg, perm_flag)

% Run nested cross-validation based on options specified in cfg
% Joseph Griffis 2024

% Set repeats to 1 if set to 0 or not set so it doesn't break
if cfg.cv.repeats == 0 || isempty(cfg.cv.repeats)
    cfg.cv.repeats = 1;
end

% Preallocate relevant fields
if perm_flag == 0
    cv_results = preallocate_cv_results(size(X,1), size(X,2), cfg);
elseif perm_flag == 1
    
    % Preallocate relevant fields
    cv_results = preallocate_cv_results_perm(size(X,1), cfg);

    % Permute relevant variables
    perm_order = randperm(length(Y));
    Y = Y(perm_order);
    if ~isempty(cfg.strat_var)
        cfg.strat_groups = cfg.strat_groups(perm_order);
    end
    if ~isempty(cfg.confounds)
        cfg.confounds = cfg.confounds(perm_order,:);
    end

end

% Run nested CV, either creating new partitions or using pre-defined partitions
if isfield(cfg.cv, 'partitions')
    cv_results = run_pre_partitioned_cv(X, Y, cfg, cv_results, perm_flag);
else
    cv_results = run_new_partition_cv(X, Y, cfg, cv_results, perm_flag);
end

% Define summary function 
if strcmp(cfg.cv.summary_type, 'mean') 
    get_summary = @nanmean;
elseif strcmp(cfg.cv.summary_type, 'median')
    get_summary = @nanmedian;
end

% Get average coefficients across outer CV loops for each repeat
if ~strcmp(cfg.model_spec, 'olsr') && perm_flag == 0
    cv_results.avg.coeff = squeeze(get_summary(cv_results.coeff, 3));
end

% Get average confound variance explained 
if ~isempty(cfg.confounds) && perm_flag == 0
    cv_results.avg.r2_confound = get_summary(cv_results.r2_confound,3);
end

% Restructure predicted and observed outcomes for outer loop test set 
if perm_flag == 0
    cv_results.all.pred_y = get_summary(cv_results.pred_y, 3); % collapse over NaNs to get predictions across all folds for each repeat
    cv_results.all.obs_y = Y; % observed scores for all patients
    if cfg.cat_Y == 1
        cv_results.all.pred_score = get_summary(cv_results.pred_score,3);
    end
end

% Get correlation (or Fisher exact test) between mean out-of-fold predictions for all patients and observed outcome, compute p-value (analogous to LESYMAP)
if cfg.cat_Y == 0 && perm_flag == 0
    [cv_results.all.corr, cv_results.all.corr_pval] = corr(get_summary(cv_results.all.pred_y, 2), Y);
    disp(['Full-sample Cross-Validation Correlation Test: r=' num2str(cv_results.all.corr) ', p=' num2str(cv_results.all.corr_pval)]);
elseif cfg.cat_Y == 1 && perm_flag == 0
    cv_results.all.confusion = confusionchart(Y, mode(cv_results.all.pred_y,2));
    [~, cv_results.all.fisher_pval, cv_results.all.fisher_stat] = fishertest(cv_results.all.confusion.NormalizedValues);
    cv_results.all.classrate(1) = cv_results.all.confusion.NormalizedValues(1,1) ./ (cv_results.all.confusion.NormalizedValues(1,2) + cv_results.all.confusion.NormalizedValues(1,1));
    cv_results.all.classrate(2) = cv_results.all.confusion.NormalizedValues(2,2) ./ (cv_results.all.confusion.NormalizedValues(2,2) + cv_results.all.confusion.NormalizedValues(2,1));    
    roc_obj = rocmetrics(cv_results.all.obs_y, get_summary(cv_results.all.pred_score,2), 1);
    cv_results.all.roc_auc = roc_obj.AUC;
    disp(['Full-sample Cross-Validation Fisher Exact Test: odds ratio=' num2str(cv_results.all.fisher_stat.OddsRatio) ', p=' num2str(cv_results.all.fisher_pval)]);
    disp(['Full-sample Cross-Validation Test Classification Accuracy (Group = -1): ' num2str(cv_results.all.classrate(1))]);
    disp(['Full-sample Cross-Validation Test Classification Accuracy (Group = +1): ' num2str(cv_results.all.classrate(2))]);
    disp(['Full-sample Cross-Validation Test Area Under ROC Curve: ' num2str(cv_results.all.roc_auc)]);
end

% Get performance measures for each repeat
if cfg.cat_Y == 0 && perm_flag == 0
    cv_results.avg.explained(:,1) = get_summary(cv_results.explained,2);
    cv_results.avg.r2_ss(:,1) = get_summary(cv_results.r2_ss,2);
    cv_results.avg.corr(:,1) = get_summary(cv_results.corr,2);
    cv_results.avg.mse(:,1) = get_summary(cv_results.mse, 2);
elseif cfg.cat_Y == 1 && perm_flag == 0
    cv_results.avg.classrate(:,1) = get_summary(cv_results.classrate(:,:,1), 2);
    cv_results.avg.classrate(:,2) = get_summary(cv_results.classrate(:,:,2), 2);
    cv_results.avg.classrate(:,3) = get_summary(cv_results.classrate(:,:,3), 2);
    cv_results.avg.roc_auc(:,1) = get_summary(cv_results.roc_auc, 2);
elseif cfg.cat_Y == 0 && perm_flag == 1
    cv_results.avg.mse = get_summary(cv_results.mse,2);
elseif cfg.cat_Y == 1 && perm_flag == 1
    cv_results.avg.roc_auc(:,1) = get_summary(cv_results.roc_auc, 2);
end

% Save CV results
if cfg.cv.save_cv_results == 1 && perm_flag == 0
    cd(cfg.out_dir);
    if isfile('cv_results.mat')
        save('cv_results.mat', 'cv_results', '-append');
    else
        save('cv_results.mat', 'cv_results', '-v7.3');
    end
end