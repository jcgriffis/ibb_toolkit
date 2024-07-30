function cv_results = run_stacked_cv(Y,cfg,cv_results,base_results,perm_flag)

% Run cross-validation using pre-defined train/test partitions
% Train/test partitions must be n_observations X n_repeats X n_folds
% Joseph Griffis 2024

for i = 1:cfg.cv.repeats % Loop over CV repeats     

    disp(['Outer CV repeats:' num2str(i) '/' num2str(cfg.cv.repeats)]);
 
    % Stack predictors across models
    X_stack = zeros(size(base_results.pred_y{1}, 1), length(base_results.pred_y));
    for k = 1:length(base_results.pred_y)
        % Predictors are predictions from each model for this iteration of repeats and folds
        if cfg.cat_Y == 0
            X_stack(:,k) = base_results.all.pred_y{k}(:,i);          
        elseif cfg.cat_Y
            if isfield(cfg, 'use_scores')
                X_stack(:,k) = base_results.all.pred_score{k}(:,i);
            elseif isfield(cfg, 'use_labels')
                X_stack(:,k) = base_results.all.pred_y{k}(:,i);
            else 
                X_stack(:,k) = base_results.all.pred_score{k}(:,i);
            end
        end
    end

    % Loop over outer folds
    for j = 1:cfg.cv.folds
        
        % Get training set
        train_set = cfg.cv.partitions.train_set(:,i,j);
        cfg.train_set = train_set;

        % Get test set
        test_set = cfg.cv.partitions.test_set(:,i,j);
        cfg.test_set = test_set;

        % Get train and test data 
        x_test = X_stack(test_set==1,:);
        y_test = Y(test_set==1);  
        x_train = X_stack(train_set==1,:);
        y_train = Y(train_set==1);

        % Apply standardization if specified
        if cfg.standardize == 1  % Convert X to z-scores
            miss_col = sum(x_train,1)==0;
            [x_train, cfg.Cx, cfg.Sx] = normalize(x_train);
            x_train(:,miss_col)=0; % Since normalize will cause columns to become NaN if they are all 0
            x_test = normalize(x_test, 'center', cfg.Cx, 'scale', cfg.Sx);
            x_test(:,miss_col)=0; % Since normalize will cause columns to become NaN if they have SD = 0 in train set
        elseif cfg.standardize == 2 % Convert Y to z-scores 
            if cfg.cat_Y == 0
                cfg.standardize = 1;
                disp('Unable to standardize categorical Y, setting standardize flag to 0');                
                [y_train, cfg.Cy, cfg.Sy] = normalize(y_train);
                y_test = normalize(y_test, 'center', cfg.Cy, 'scale', cfg.Sy);
            end
        elseif cfg.standardize == 3 % Convert X and Y to z-scores
            miss_col = sum(x_train,1)==0;
            [x_train, cfg.Cx, cfg.Sx] = normalize(x_train);
            x_train(:,miss_col)=0; % Since normalize will cause columns to become NaN if they are all 0
            x_test = normalize(x_test, 'center', cfg.Cx, 'scale', cfg.Sx);
            x_test(:,miss_col)=0; % Since normalize will cause columns to become NaN if they have SD = 0 in train set
            if cfg.cat_Y == 0
                disp('Unable to standardize categorical Y, setting standardize flag to 1');
                cfg.standardize = 1;
                [y_train, cfg.Cy, cfg.Sy] = normalize(y_train);
                y_test = normalize(y_test, 'center', cfg.Cy, 'scale', cfg.Sy);      
            end
        end
    
        % Run inner loop and get cross-validation results
        switch cfg.model_spec 
            case 'plsr'
                cv_results = run_plsr_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag);
                clear y_test x_test y_test x_train y_train train_set test_set
            case 'pls_da'
                cv_results = run_plsda_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag);
                clear y_test x_test y_test x_train y_train train_set test_set
            case {'ridge', 'lasso', 'rlinsvr'} 
                cv_results = run_reg_linear_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag);
                clear y_test x_test y_test x_train y_train train_set test_set
            case {'logistic_ridge', 'logistic_lasso', 'rlinsvc'}
                cv_results = run_class_linear_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag);
                clear y_test x_test y_test x_train y_train train_set test_set        
            case {'kernsvr', 'linsvr'}
                cv_results = run_svr_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag);
                clear y_test x_test y_test x_train y_train train_set test_set
            case {'kernsvc', 'linsvc'}
                cv_results = run_svc_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag);
                clear y_test x_test y_test x_train y_train train_set test_set      
            case 'censemble'
                cv_results = run_censemble_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag);
                clear y_test x_test y_test x_train y_train train_set test_set                      
            case 'rensemble'
                cv_results = run_rensemble_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag);   
                clear y_test x_test y_test x_train y_train train_set test_set                      
            case 'olsr'
                cv_results = run_ols_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag);
                clear y_test x_test y_test x_train y_train train_set test_set     
            case 'mean'
                cv_results = get_mean_prediction(x_test, y_test, cfg, cv_results, i, j);
                clear y_test x_test y_test x_train y_train train_set test_set           
            case 'mean_score'
                cv_results = get_mean_score_prediction(x_train, y_train, x_test, y_test, cfg, cv_results, i, j);
        end
    end    
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

% Get predictor matrix
cv_results.X_stack = X_stack;

% Save CV results
if cfg.cv.save_cv_results == 1 && perm_flag == 0
    cd(cfg.out_dir);
    if isfile('cv_results.mat')
        save('cv_results.mat', 'cv_results', 'cfg', '-append');
    else
        save('cv_results.mat', 'cv_results', 'cfg', '-v7.3');
    end
end

end
