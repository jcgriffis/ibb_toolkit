function [model_results] = run_svc_modeling(X, Y, cfg, model_results)

% Use K-fold cross-validation to tune SVC hyper-parameters and fit model
% Joseph Griffis 2023

% Get costs
S.ClassNames = [1, -1];
S.ClassificationCosts = cfg.cost;

if cfg.optimize_hyperparams == 1

    % Get hyper-parameter options and cross-validation partitions
    model_results.hp_opt = get_hp_opts(cfg, length(Y));
    
    % Run optimization
    alpha = []; % sometimes model returns empty results despite convergence
    while isempty(alpha)
        mdl_final = fitcsvm(X,Y,'KernelFunction', cfg.kernel,...
            'OptimizeHyperparameters', {'BoxConstraint', 'KernelScale'},...
            'HyperparameterOptimizationOptions', model_results.hp_opt, ...
            'Cost', S);
        if isnan(mdl_final.Bias)
            alpha = [];
        else
            alpha = mdl_final.Alpha;
        end
    end
else

    alpha = []; % sometimes model returns empty results despite convergence
    while isempty(alpha)
        mdl_final = fitcsvm(X,Y,'KernelFunction', cfg.kernel,...
            'Cost', S);
        if isnan(mdl_final.Bias)
            alpha = [];
        else
            alpha = mdl_final.Alpha;
        end
    end

end

% Get optimal hyper-parameters
model_results.C = mdl_final.BoxConstraints(1);
model_results.gamma = mdl_final.KernelParameters.Scale;
model_results.alpha = mdl_final.Alpha;
clear mdl 

% Compute betas using sensitivity mapping method (Zhang et al., 2014 - Human Brain Mapping)
if strcmp(cfg.model_spec, 'kernsvc')
    model_results.coeff(:,1) = (mdl_final.Alpha.'*(mdl_final.SupportVectors.*mdl_final.SupportVectorLabels));
else
    model_results.coeff(:,1) = mdl_final.Beta;
end

% Get final model classification rates
[pred_y, pred_score] = predict(mdl_final, X);

% Get ROC AUC 
roc_obj = rocmetrics(Y, pred_score(:,2), 1);    
model_results.roc_auc = roc_obj.AUC;

% Get optimal threshold
opt_thresh = get_optimal_threshold(Y, roc_obj, cfg.cost);

% Assign labels using optimized threshold
pred_y(:) = 0;
pred_y(pred_score(:,2) >= opt_thresh) = 1;
pred_y(pred_score(:,2) < opt_thresh) = -1;

% Get final model classification rates
model_results.classrate(2) = numel(intersect(find(pred_y==1), find(Y==1)))./numel(find(Y==1));
model_results.classrate(1) = numel(intersect(find(pred_y==-1), find(Y==-1)))./numel(find(Y==-1));
model_results.classrate(3) = numel(find(pred_y==Y))./numel(find(Y));

% Get ROC AUC 
roc_obj = rocmetrics(Y, pred_score(:,2), 1);    
model_results.roc_auc = roc_obj.AUC;    

% Save predicted and observed
model_results.pred_y = pred_y; % Predicted outcomes
model_results.pred_score = pred_score(:,2); % Prediction scores
model_results.obs_y = Y; % Observed outcomes
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