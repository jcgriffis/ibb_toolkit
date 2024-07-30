function [model_results] = run_class_linear_modeling(X, Y, cfg, model_results)

% Use K-fold cross-validation to tune regularized regression hyper-parameters
% and fit final model to full dataset. 
% Joseph Griffis 2023

% Get costs
S.ClassNames = [1, -1];
S.ClassificationCosts = cfg.cost;

if cfg.optimize_hyperparams == 1
    
    %%%% Do hyperparameter optimization and fit model
    
    % CV hyper-parameter optimization and model fitting
    if cfg.hp_opt.bayes_opt == 0
        
        % Define lambda range and preallocate output
        lambda_range = logspace(-5, 5, 30) ./ size(X,1);
        lambdas = zeros(cfg.hp_opt.repeats,1);
        
        % Remove dot indexing for parallel processing
        reg_type = cfg.reg_type;
        learner = cfg.learner;
        repeats = cfg.hp_opt.repeats;
        cv_part = get_hp_opts(cfg, length(Y));

        if cfg.parallel == 1
            parfor i = 1:cfg.hp_opt.repeats
                
                % Repartition
                my_part = repartition(cv_part);
                
                % Get hyper-parameter options and cross-validation partitions
                disp(['Hyper-parameter optimization iteration:' num2str(i) '/' num2str(repeats)]) 
                
                mdl = fitclinear(X,Y,...
                'ObservationsIn', 'rows',...
                'Regularization', reg_type,...
                'Learner', learner,...
                'Lambda', lambda_range,...
                'CVPartition', my_part,...
                'Cost', S,...
                'BetaTolerance', 0);
                            
                % Get optimal lambda for each iteration
                loss = kfoldLoss(mdl);
                min_idx = lambda_range(loss == min(loss));
                lambdas(i) = min_idx(1);
            end
        else
            for i = 1:cfg.hp_opt.repeats
        
                disp(['Hyper-parameter optimization iteration:' num2str(i) '/' num2str(cfg.hp_opt.repeats)]) 
                
                % Repartition train/test
                cv_part = repartition(cv_part);

                mdl = fitclinear(X,Y,...
                'ObservationsIn', 'rows',...
                'Regularization', reg_type,...
                'Learner', learner,...
                'Lambda', lambda_range,...
                'CVPartition', cv_part,...
                'Cost', S,...
                'BetaTolerance', 0);
                            
                % Get optimal lambda for each iteration
                loss = kfoldLoss(mdl);
                min_idx = lambda_range(loss == min(loss));
                lambdas(i) = min_idx(1);
                
            end
        end

        % Choose optimal lambda
        opt_l = mean(lambdas);
        
        disp(['Optimal lambda = ' num2str(opt_l) ', fitting model with optimal lambda']);
    
        % Fit model 
        mdl_final = fitclinear(X,Y,...
        'ObservationsIn', 'rows',...
        'Regularization', cfg.reg_type,...
        'Lambda', opt_l,...
        'Learner', cfg.learner,...
        'Cost', S,...
        'BetaTolerance', 0);
           
    else
            disp('Optimizing lambda with Bayesian optmization');
            
            model_results.hp_opt = get_hp_opts(cfg, length(Y));
     
            mdl_final = fitclinear(X,Y,...
            'ObservationsIn', 'rows',...
            'Regularization', cfg.reg_type,...
            'OptimizeHyperparameters', {'Lambda'},...
            'Learner', cfg.learner,...
            'HyperparameterOptimizationOptions', model_results.hp_opt,...
            'Cost', S,...
            'BetaTolerance', 0);
    
    end
else
    mdl_final = fitrlinear(X,Y,...
    'ObservationsIn', 'rows',...
    'Regularization', cfg.reg_type,...
    'Learner', cfg.learner,...
    'Cost', S,...
    'BetaTolerance', 0);     
end

% Get final model classification rates
[pred_y, pred_score] = predict(mdl_final, X, 'ObservationsIn', 'rows');

% Get ROC AUC 
roc_obj = rocmetrics(Y, pred_score(:,2), 1);    
model_results.roc_auc = roc_obj.AUC;

% Get optimal threshold
opt_thresh = get_optimal_threshold(Y, roc_obj, cfg.cost);

% Assign labels using optimized threshold
pred_y(:) = 0;
pred_y(pred_score(:,2) >= opt_thresh) = 1;
pred_y(pred_score(:,2) < opt_thresh) = -1;

% Get classification rates
model_results.classrate(2) = numel(intersect(find(pred_y==1), find(Y==1)))./numel(find(Y==1));
model_results.classrate(1) = numel(intersect(find(pred_y==-1), find(Y==-1)))./numel(find(Y==-1));
model_results.classrate(3) = numel(find(pred_y==Y))./numel(find(Y));

% Get beta weights
model_results.coeff = mdl_final.Beta;
model_results.beta0 = mdl_final.Bias;
model_results.Lambda = mdl_final.Lambda;

% Store predicted and observed Y
model_results.obs_y = Y;
model_results.pred_y = pred_y;
model_results.pred_score = pred_score(:,2);
model_results.opt_score_thresh = opt_thresh;

% Save final model if indicated
if cfg.save_model_results == 1
    if isfile('model_results.mat')
       save('model_results.mat', 'mdl_final', '-append');
    else
       save('model_results.mat', 'mdl_final', '-v7.3');
    end
end

end