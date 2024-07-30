function [cv_results] = run_censemble_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag)

% Use K-fold cross-validation to tune classification ensemble hyper-parameters and fit model in nested CV scheme
% Stein Acker 2024
% (heavily copied from run_svc_modeling_cv.m by Joseph Griffis)

% Get costs
S.ClassNames = [1, -1];
S.ClassificationCosts = cfg.cost;

if cfg.optimize_hyperparams == 1

    % Get hyper-parameter options and cross-validation partitions
    model_results.hp_opt = get_hp_opts(cfg, length(y_train));

    % Run optimization
    if iscellstr(cfg.hp_opt.to_optimize)
        if contains('Method', cfg.hp_opt.to_optimize)
            mdl = fitcensemble(x_train,y_train,'Learners', 'tree',...
                'OptimizeHyperparameters', cfg.hp_opt.to_optimize,...
                'HyperparameterOptimizationOptions', model_results.hp_opt,...
                'Cost', S);

        else
            mdl = fitcensemble(x_train,y_train,'Learners', 'tree',...
                'Method', cfg.method,...
                'Options', statset('UseParallel',cfg.parallel),...
                'Cost', S);
        end
    elseif isa(cfg.hp_opt.to_optimize(1), 'optimizableVariable')
        if cfg.hp_opt.to_optimize(1).Optimize == 1
            mdl = fitcensemble(x_train,y_train,'Learners', 'tree',...
                'OptimizeHyperparameters', cfg.hp_opt.to_optimize,...
                'HyperparameterOptimizationOptions', model_results.hp_opt,...
                'Cost', S);
        else
            mdl = fitcensemble(x_train,y_train,'Learners', 'tree',...
                'Method',cfg.method,...
                'OptimizeHyperparameters', cfg.hp_opt.to_optimize,...
                'HyperparameterOptimizationOptions', model_results.hp_opt,...
                'Cost', S);
        end
    end
else
    mdl = fitcensemble(x_train,y_train,'Learners', 'tree',...
        'Method', cfg.method,...
        'Options', statset('UseParallel',cfg.parallel),...
        'Cost', S);
end

% Get optimal threshold from training set
[~, pred_score] = predict(mdl, x_train);
roc_obj = rocmetrics(y_train, pred_score(:,2), 1);    
opt_thresh = get_optimal_threshold(y_train, roc_obj, cfg.cost);

% Get predictions and coefficients
[pred_y, pred_score] = predict(mdl, x_test);

% Assign labels using optimized threshold
pred_y(:) = 0;
pred_y(pred_score(:,2) >= opt_thresh) = 1;
pred_y(pred_score(:,2) < opt_thresh) = -1;

% Store predictions and observations
cv_results.pred_y(cfg.test_set==1,i,j) = pred_y; % Predicted outcome for test set
cv_results.obs_y(cfg.test_set==1,i,j) = y_test; % Observed outcome for test set
roc_obj = rocmetrics(y_test, pred_score(:,2), 1);    
cv_results.roc_auc(i,j) = roc_obj.AUC;   

if perm_flag == 0

    % Get predictor importances
    if strcmp(mdl.ModelParameters.Method, 'Bag')
        cv_results.coeff(:,i,j) = oobPermutedPredictorImportance(mdl, 'Options', statset('UseParallel', cfg.parallel));
    else
        cv_results.coeff(:,i,j) = predictorImportance(mdl);
    end

    % Get method
    cv_results.method{i,j} = mdl.ModelParameters.Method;
    cv_results.n_learn(i,j) = mdl.ModelParameters.NLearn;

    % Get scores and performance measures
    cv_results.pred_score(cfg.test_set==1,i,j) = pred_score(:,2);
    cv_results.classrate(i,j,3) = numel(find(pred_y==y_test))./numel(find(y_test));
    cv_results.classrate(i,j,2) = numel(intersect(find(pred_y==1), find(y_test==1)))./numel(find(y_test==1));
    cv_results.classrate(i,j,1) = numel(intersect(find(pred_y==-1), find(y_test==-1)))./numel(find(y_test==-1));
end

end