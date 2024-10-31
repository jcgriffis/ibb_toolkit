function [cv_results] = run_class_linear_modeling_cv(x_train, y_train, x_test, y_test, i, j, cfg, cv_results, perm_flag)

% Use K-fold cross-validation to tune regularized regression hyper-parameters
% and predict scores for test set
% Joseph Griffis 2023

% Get classification costs
S.ClassNames = [1, -1];
S.ClassificationCosts = cfg.cost;

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
                
                % Repartition
                my_part = repartition(cv_part);
                
                % Get hyper-parameter options and cross-validation partitions
                disp(['Hyper-parameter optimization iteration:' num2str(k) '/' num2str(repeats)]) 
                
                mdl = fitclinear(x_train,y_train,...
                'ObservationsIn', 'rows',...
                'Regularization', reg_type,...
                'Learner', learner,...
                'Lambda', lambda_range,...
                'CVPartition', my_part,...
                'BetaTolerance', 0,...
                'Cost', S);
                            
                % Get optimal lambda for each iteration
                loss = kfoldLoss(mdl);
                min_idx = lambda_range(loss == min(loss));
                lambdas(k) = min_idx(1);
            end
        else
            for k = 1:cfg.hp_opt.repeats
                
                my_part = repartition(cv_part);

                disp(['Hyper-parameter optimization iteration:' num2str(k) '/' num2str(cfg.hp_opt.repeats)]) 
                
                mdl = fitclinear(x_train,y_train,...
                'ObservationsIn', 'rows',...
                'Regularization', reg_type,...
                'Learner', learner,...
                'Lambda', lambda_range,...
                'CVPartition', my_part,...
                'BetaTolerance', 0, ...
                'Cost', S);
                            
                % Get optimal lambda for each iteration
                loss = kfoldLoss(mdl);
                lambdas(k) = lambda_range(loss == min(loss));
                
            end
        end

        % Choose optimal lambda
        opt_l = mean(lambdas);
        
        disp(['Optimal lambda = ' num2str(opt_l) ', fitting model with optimal lambda']);
    
        % Fit model 
        mdl = fitclinear(x_train,y_train,...
        'ObservationsIn', 'rows',...
        'Regularization', cfg.reg_type,...
        'Lambda', opt_l,...
        'Learner', cfg.learner,...
        'BetaTolerance', 0,...
        'Cost', S);
           
    else
            disp('Optimizing lambda with Bayesian optmization');
           
            cv_results.hp_opt = get_hp_opts(cfg, length(y_train));
      
            mdl = fitclinear(x_train,y_train,...
            'ObservationsIn', 'rows',...
            'Regularization', cfg.reg_type,...
            'OptimizeHyperparameters', {'Lambda'},...
            'Learner', cfg.learner,...
            'HyperparameterOptimizationOptions', cv_results.hp_opt,...
            'BetaTolerance', 0,...
            'Cost',S);
    
    end
else
    mdl = fitclinear(x_train,y_train,...
    'ObservationsIn', 'rows',...
    'Regularization', cfg.reg_type,...
    'Learner', cfg.learner,...
    'BetaTolerance', 0,...
    'Cost', S);     
    opt_l = mdl.Lambda;
end

% Get optimal threshold from training set
[~, pred_score] = predict(mdl, x_train, 'ObservationsIn', 'rows');
roc_obj = rocmetrics(y_train, pred_score(:,2), 1);    
opt_thresh = get_optimal_threshold(y_train, roc_obj, cfg.cost);

% Get predictions and coefficients
[pred_y, pred_score] = predict(mdl, x_test, 'ObservationsIn', 'rows');

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
    cv_results.pred_score(cfg.test_set==1,i,j) = pred_score(:,2);    
    cv_results.classrate(i,j,2) = numel(intersect(find(pred_y==1), find(y_test==1)))./numel(find(y_test==1));
    cv_results.classrate(i,j,1) = numel(intersect(find(pred_y==-1), find(y_test==-1)))./numel(find(y_test==-1));    
    cv_results.classrate(i,j,3) = numel(find(pred_y==y_test))./numel(find(y_test));
    cv_results.coeff(:,i,j) = mdl.Beta; % beta coefficients    
    cv_results.beta_0(:,i,j) = mdl.Bias; % intercept
    cv_results.lambda(i,j) = opt_l; % optimal lambda from training set
    cv_results.opt_thresh(i,j) = opt_thresh;
end

end
