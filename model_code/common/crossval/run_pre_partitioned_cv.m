function cv_results = run_pre_partitioned_cv(X,Y,cfg,cv_results,perm_flag)

% Run cross-validation using pre-defined train/test partitions
% Train/test partitions must be n_observations X n_repeats X n_folds
% Joseph Griffis 2024

for i = 1:cfg.cv.repeats % Loop over CV repeats     

    disp(['Outer CV repeats:' num2str(i) '/' num2str(cfg.cv.repeats)]);
    
    % Loop over outer folds
    for j = 1:cfg.cv.folds
        
        % Get training set
        train_set = cfg.cv.partitions.train_set(:,i,j);
        cfg.train_set = train_set;
        if perm_flag == 0
            cv_results.train_set(:,i,j) = train_set;
        end
        
        % Get test set
        test_set = cfg.cv.partitions.test_set(:,i,j);
        cfg.test_set = test_set;
        if perm_flag == 0
            cv_results.test_set(:,i,j) = test_set;
        end

        % Get train and test data 
        x_test = X(test_set==1,:);
        y_test = Y(test_set==1);  
        x_train = X(train_set==1,:);
        y_train = Y(train_set==1);
        
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
        if cfg.standardize == 1  % Convert X to z-scores
            miss_col = sum(x_train,1)==0;
            [x_train, cfg.Cx, cfg.Sx] = normalize(x_train);
            x_train(:,miss_col)=0; % Since normalize will cause columns to become NaN if they are all 0
            x_test = normalize(x_test, 'center', cfg.Cx, 'scale', cfg.Sx);
            x_test(:,miss_col)=0; % Since normalize will cause columns to become NaN if they have SD = 0 in train set
        elseif cfg.standardize == 2 % Convert Y to z-scores 
            if cfg.cat_Y == 1
                cfg.standardize = 0;
                disp('Unable to standardize categorical Y, setting standardize flag to 0');    
            else
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
                cfg.standardize = 0;
            else
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
                cv_results = run_censemble_modeling(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag);
                clear y_test x_test y_test x_train y_train train_set test_set                      
            case 'rensemble'
                cv_results = run_rensemble_modeling(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag);   
                clear y_test x_test y_test x_train y_train train_set test_set                      
            case 'olsr'
                cv_results = run_ols_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag);
                clear y_test x_test y_test x_train y_train train_set test_set              
        end
        if isfield(cfg, 'Cx')
            cfg = rmfield(cfg, {'Cx', 'Cy'});
        end
        if isfield (cfg, 'Cy')
            cfg = rmfield(cfg, {'Cx', 'Cy'});
        end
    end    
end
