function cv_results = run_new_partition_cv(X, Y, cfg, cv_results, perm_flag)

% Run cross-validation without pre-defined train/test partitions
% Joseph Griffis 2024

% Pre-generate all partitions if stacking is desired (i.e., to avoid potential changes in RNG state)
if cfg.cv.stacked_model == 1 && strcmp(cfg.cv.cv_type, 'KFold') 
        all_partitions = cell(cfg.cv.repeats,1);    
    for i = 1:cfg.cv.repeats
        if isempty(cfg.strat_var)
            all_partitions{i} = cvpartition(length(Y), 'KFold', cfg.cv.folds);
        else
            all_partitions{i} = cvpartition(cfg.strat_groups,'KFold', cfg.cv.folds, 'Stratify', true);
        end
    end
end

for i = 1:cfg.cv.repeats % Loop over CV repeats     

    disp(['Outer CV repeats:' num2str(i) '/' num2str(cfg.cv.repeats)]);

    % Get partition of outer train/test folds for this repeat
    if cfg.cv.stacked_model == 1 && strcmp(cfg.cv.cv_type, 'KFold')
        outer_folds = all_partitions{i};
    else % Generate outer folds on first iteration if stacking not desired
        if i == 1
            if ~isempty(cfg.strat_var) % Get groups based on stratification variable and define folds
                outer_folds = cvpartition(cfg.strat_groups,'KFold', cfg.cv.folds, 'Stratify', true);
            else % Define folds
                switch cfg.cv.cv_type
                    case 'KFold'
                        outer_folds = cvpartition(length(Y),'KFold', cfg.cv.folds);
                    case 'LOOCV'
                        outer_folds = cvpartition(length(Y), 'LeaveOut');
                end
            end
        else % repartition on later iterations
            outer_folds = outer_folds.repartition;
        end
    end

    % Loop over outer folds
    for j = 1:outer_folds.NumTestSets   
        
        % Get training set
        train_set = training(outer_folds, j);
        cfg.train_set = train_set;
        if perm_flag == 0
            cv_results.train_set(:,i,j) = train_set;
        end
        
        % Get test set
        test_set = test(outer_folds, j);
        cfg.test_set = test_set;
        if perm_flag == 0
            cv_results.test_set(:,i,j) = test_set;
        end

        % Get train and test data 
        x_test = X(test_set,:);
        y_test = Y(test_set);  
        x_train = X(train_set,:);
        y_train = Y(train_set);
        
        % Do confound regression if specified
        if cfg.cat_Y == 0 % Only do this for continuous outcomes
            if ~isempty(cfg.confounds) % Confound regression                
                % Regress Y on confounds and get residual Y
                disp('Regressing confounds out of Y...');
                [y_train, cv_results.r2_confound(i,j), cv_results.coeff_confound(:,i,j)] = regress_confounds_from_Y(y_train, cfg.confounds(train_set,:));
                disp(['Confound model explained ' num2str(round(cv_results.r2_confound(i,j), 4)*100) '% variance in Y.']);
                disp('Continuing analysis with Y residuals as target variable');    
                % Fit training confound model to test set and get residuals (new Y)
                y_pred_confound = [ones(length(y_test), 1) cfg.confounds(test_set,:)]*cv_results.coeff_confound(:,i,j);
                y_test = y_test - y_pred_confound;
                clear y_pred_confound
            end
        end

        % Apply standardization if specified
        if cfg.standardize > 0  
           [cfg, x_train, x_test, y_train, y_test] = rescale_dataset(cfg, x_train, x_test, y_train, y_test);
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
        end
    end    
end
end
