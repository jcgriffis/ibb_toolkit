function [cv_results] = run_rensemble_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag)

% Use K-fold cross-validation to tune classification ensemble hyper-parameters and fit model in nested CV scheme
% Stein Acker 2024
% (heavily copied from run_svc_modeling_cv.m by Joseph Griffis, modified by Joseph Griffis)

if cfg.optimize_hyperparams == 1

    % Get hyper-parameter options and cross-validation partitions
    model_results.hp_opt = get_hp_opts(cfg, length(y_train));

    % Run optimization
    if iscellstr(cfg.hp_opt.to_optimize)
        if contains('Method', cfg.hp_opt.to_optimize)
            mdl = fitrensemble(x_train,y_train,'Learners', 'tree',...
                'OptimizeHyperparameters', cfg.hp_opt.to_optimize,...
                'HyperparameterOptimizationOptions', model_results.hp_opt);

        else
            mdl = fitrensemble(x_train,y_train,'Learners', 'tree',...
                'Method', cfg.method,...
                'Options', statset('UseParallel',cfg.parallel));
        end
    elseif isa(cfg.hp_opt.to_optimize(1), 'optimizableVariable')
        if cfg.hp_opt.to_optimize(1).Optimize == 1
            mdl = fitrensemble(x_train,y_train,'Learners', 'tree',...
                'OptimizeHyperparameters', cfg.hp_opt.to_optimize,...
                'HyperparameterOptimizationOptions', model_results.hp_opt);
        else
            mdl = fitrensemble(x_train,y_train,'Learners', 'tree',...
                'Method',cfg.method,...
                'OptimizeHyperparameters', cfg.hp_opt.to_optimize,...
                'HyperparameterOptimizationOptions', model_results.hp_opt);
        end
    end
else
    mdl = fitrensemble(x_train,y_train,'Learners', 'tree',...
        'Method', cfg.method,...
        'Options', statset('UseParallel',cfg.parallel));
end

% Predict test cases
y_pred = predict(mdl, x_test, 'ObservationsIn', 'rows');

% Put predictions and observations back in original units if needed
if isfield(cfg, 'Cy')
    y_pred = (y_pred .* cfg.Sy) + cfg.Cy;
    y_test = (y_test .* cfg.Sy) + cfg.Cy;
    y_train = (y_train .* cfg.Sy) + cfg.Cy;
end

% Store predictions and observations
cv_results.pred_y(cfg.test_set==1,i,j) = y_pred; % Predicted outcome for test set
cv_results.obs_y(cfg.test_set==1,i,j) = y_test; % Observed outcome for test set
cv_results.mse(i,j) = get_mse(y_test, y_pred); % MSE

% Store other model information
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

    % Get model performance estimates
    cv_results.explained(i,j) = get_explained_variance(y_test, y_pred); % explained variance
    cv_results.r2_ss(i,j) = get_model_r2_train(y_test, y_pred, y_train); % sum-of-squares R-squared
    cv_results.corr(i,j) = corr(y_test, y_pred); % Correlation of predicted and observed    
end

end
