function [model_results] = run_censemble_modeling(X, Y, cfg, model_results)

% Use K-fold cross-validation to tune classification ensemble hyper-parameters and fit model
% Stein Acker 2024
% (heavily copied from run_svc_modeling.m by Joseph Griffis; modified by Joseph Griffis)

% Get costs
S.ClassNames = [1, -1];
S.ClassificationCosts = cfg.cost;

if cfg.optimize_hyperparams == 1

    % Get hyper-parameter options and cross-validation partitions
    model_results.hp_opt = get_hp_opts(cfg, length(Y));

    % Run optimization
    if iscellstr(cfg.hp_opt.to_optimize)
        if contains('Method', cfg.hp_opt.to_optimize)
            mdl_final = fitcensemble(X,Y,'Learners', 'tree',...
                'OptimizeHyperparameters', cfg.hp_opt.to_optimize,...
                'HyperparameterOptimizationOptions', model_results.hp_opt,...
                'Cost', S);

        else
            mdl_final = fitcensemble(X,Y,'Learners', 'tree',...
                'Method', cfg.method,...
                'Options', statset('UseParallel',cfg.parallel),...
                'Cost', S);
        end
    elseif isa(cfg.hp_opt.to_optimize(1), 'optimizableVariable')
        if cfg.hp_opt.to_optimize(1).Optimize == 1
            mdl_final = fitcensemble(X,Y,'Learners', 'tree',...
                'OptimizeHyperparameters', cfg.hp_opt.to_optimize,...
                'HyperparameterOptimizationOptions', model_results.hp_opt,...
                'Cost', S);
        else
            mdl_final = fitcensemble(X,Y,'Learners', 'tree',...
                'Method',cfg.method,...
                'OptimizeHyperparameters', cfg.hp_opt.to_optimize,...
                'HyperparameterOptimizationOptions', model_results.hp_opt, ...
                'Cost', S);
        end
    end
else
    mdl_final = fitcensemble(X,Y,'Learners', 'tree',...
        'Method', cfg.method,...
        'Options', statset('UseParallel',cfg.parallel),...
        'Cost', S);
end

% Get optimal hyper-parameters
model_results.method = mdl_final.ModelParameters.Method;
model_results.n_learn = mdl_final.ModelParameters.NLearn;

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

% Save predicted and observed
model_results.pred_y = pred_y; % Predicted outcomes
model_results.pred_score = pred_score(:,2); % Prediction scores
model_results.obs_y = Y; % Observed outcomes
model_results.opt_score_thresh = opt_thresh;

% Get predictor importances
if strcmp(mdl_final.ModelParameters.Method, 'Bag')
    model_results.coeff = oobPermutedPredictorImportance(mdl_final, 'Options', statset('UseParallel', cfg.parallel));
else
    model_results.coeff = predictorImportance(mdl_final);
end

% Save final model if indicated
if cfg.save_model_results == 1
    if isfile('model_results.mat')
       save('model_results.mat', 'mdl_final', '-append');
    else
       save('model_results.mat', 'mdl_final', '-v7.3');
    end
end

end