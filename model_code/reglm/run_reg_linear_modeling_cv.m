function [cv_results] = run_reg_linear_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag)

% Use K-fold cross-validation to tune regularized regression hyper-parameters
% and predict scores for test set
% Joseph Griffis 2023

if cfg.optimize_hyperparams == 1
    
    %%%% Do hyperparameter optimization and fit model
    
    % CV hyper-parameter optimization and model fitting
    if cfg.hp_opt.bayes_opt == 0
        
        % Define lambda range and preallocate output
        lambda_range = logspace(-5, 5, 30) ./ size(x_train,1);
        lambdas = zeros(cfg.hp_opt.repeats,1);
        
        % Remove dot indexing for parallel processing
        reg_type = cfg.reg_type;
        learner = cfg.learner;
        repeats = cfg.hp_opt.repeats;
        cv_part = get_hp_opts(cfg, length(y_train));

        if cfg.parallel == 1
            parfor k = 1:cfg.hp_opt.repeats
    
                % Pre-define loss as NaNs, occasionally kFoldLoss returns NaN values and breaks the routine
                loss = ones(length(lambda_range),1).*NaN; 
                
                while sum(isnan(loss)) > 0
                    % Repartition
                    my_part = repartition(cv_part);
                    
                    % Get hyper-parameter options and cross-validation partitions
                    disp(['Hyper-parameter optimization iteration:' num2str(k) '/' num2str(repeats)]) 
                    
                    mdl = fitrlinear(x_train,y_train,...
                    'ObservationsIn', 'rows',...
                    'Regularization', reg_type,...
                    'Learner', learner,...
                    'Lambda', lambda_range,...
                    'CVPartition', my_part,...
                    'BetaTolerance', 0);
                                
                    % Get optimal lambda for each iteration
                    loss = kfoldLoss(mdl);
                end
                min_idx = lambda_range(loss == min(loss));
                lambdas(k) = min_idx(1);
            end
        else
            for k = 1:cfg.hp_opt.repeats

                % Pre-define loss as NaNs, occasionally kFoldLoss returns NaN values and breaks the routine
                loss = ones(length(lambda_range),1).*NaN; 
                
                while sum(isnan(loss)) > 0
                    disp(['Hyper-parameter optimization iteration:' num2str(k) '/' num2str(cfg.hp_opt.repeats)]) 
                    
                    % Repartition train/test
                    cv_part = repartition(cv_part);
    
                    mdl = fitrlinear(x_train,y_train,...
                    'ObservationsIn', 'rows',...
                    'Regularization', reg_type,...
                    'Learner', learner,...
                    'Lambda', lambda_range,...
                    'CVPartition', cv_part,...
                    'BetaTolerance', 0);
                                
                    % Get optimal lambda for each iteration
                    loss = kfoldLoss(mdl);
                end
                min_idx = lambda_range(loss == min(loss));
                lambdas(k) = min_idx(1);
                
            end
        end

        % Choose optimal lambda
        opt_l = mean(lambdas);
        
        disp(['Optimal lambda = ' num2str(opt_l) ', fitting model with optimal lambda']);
    
        % Fit model 
        mdl = fitrlinear(x_train,y_train,...
        'ObservationsIn', 'rows',...
        'Regularization', cfg.reg_type,...
        'Lambda', opt_l,...
        'Learner', cfg.learner,...
        'BetaTolerance', 0);
           
    else
            disp('Optimizing lambda with Bayesian optmization');

            cv_results.hp_opt = get_hp_opts(cfg, length(y_train));

            mdl = fitrlinear(x_train,y_train,...
            'ObservationsIn', 'rows',...
            'Regularization', cfg.reg_type,...
            'OptimizeHyperparameters', {'Lambda'},...
            'Learner', cfg.learner,...
            'HyperparameterOptimizationOptions', cv_results.hp_opt,...
            'BetaTolerance', 0);
    
    end
else
    if ~isfield(cfg, 'lambda')
        mdl = fitrlinear(x_train,y_train,...
        'ObservationsIn', 'rows',...
        'Regularization', cfg.reg_type,...
        'Learner', cfg.learner,...
        'BetaTolerance', 0);
    else
        mdl = fitrlinear(x_train,y_train,...
        'ObservationsIn', 'rows',...
        'Regularization', cfg.reg_type,...
        'Learner', cfg.learner,...
        'Lambda', cfg.lambda,...
        'BetaTolerance', 0);    
    end
end

% Predict test cases
y_pred = predict(mdl, x_test, 'ObservationsIn', 'rows');

% Put predictions and observations back in original units if needed
if isfield(cfg, 'Cy')
    y_pred = (y_pred .* cfg.Sy) + cfg.Cy;
    y_test = (y_test .* cfg.Sy) + cfg.Cy;
end

% Store predictions and observations
cv_results.pred_y(cfg.test_set==1,i,j) = y_pred; % Predicted outcome for test set
cv_results.obs_y(cfg.test_set==1,i,j) = y_test; % Observed outcome for test set
cv_results.mse(i,j) = get_mse(y_test, y_pred); % MSE

% Store other model information
if perm_flag == 0
    cv_results.coeff(:,i,j) = mdl.Beta; % beta coefficients    
    cv_results.explained(i,j) = get_explained_variance(y_test, y_pred); % explained variance
    cv_results.r2_ss(i,j) = get_model_r2(y_test, y_pred); % sum-of-squares R-squared
    cv_results.corr(i,j) = corr(y_test, y_pred); % Correlation of predicted and observed    
    cv_results.beta_0(i,j) = mdl.Bias; % intercept
    cv_results.lambda(i,j) = opt_l; % optimal lambda from training set
end

end