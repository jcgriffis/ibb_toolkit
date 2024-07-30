function [model_results] = run_svr_modeling(X, Y, cfg, model_results)

% Use K-fold cross-validation to tune SVR hyper-parameters
% Joseph Griffis 2023

% Optimize hyper-parameters if indicated
if cfg.optimize_hyperparams == 1

    % Get hyper-parameter options and cross-validation partitions
    model_results.hp_opt = get_hp_opts(cfg, length(Y));
    
    % Run optimization
    alpha = []; % sometimes model returns empty results despite convergence
    while isempty(alpha)
        mdl_final = fitrsvm(X,Y,'KernelFunction', cfg.kernel,...
            'OptimizeHyperparameters', {'BoxConstraint', 'KernelScale'},...
            'HyperparameterOptimizationOptions', model_results.hp_opt);
        if isnan(mdl_final.Bias)
            alpha = [];
        else
            alpha = mdl_final.Alpha;
        end
    end
else

    alpha = []; % sometimes model returns empty results despite convergence
    while isempty(alpha)
        mdl_final = fitrsvm(X,Y,'KernelFunction', cfg.kernel);
        if isnan(mdl_final.Bias)
            alpha = [];
        else
            alpha = mdl_final.Alpha;
        end
    end

end

% Get optimal hyper-parameters
model_results.C = mdl_final.BoxConstraints(1);
model_results.epsilon = mdl_final.Epsilon;
model_results.gamma = mdl_final.KernelParameters.Scale;

% Compute betas using sensitivity mapping method (Zhang et al., 2014 - Human Brain Mapping)
if strcmp(cfg.model_spec, 'kernsvr')
    model_results.coeff(:,1) = (mdl_final.Alpha.'*mdl_final.SupportVectors);
else
    model_results.coeff(:,1) = mdl_final.Beta;
end
% Get final model R-squared
model_results.pred_y = predict(mdl_final, X); % Predicted outcome for test set

% Put predictions and observations back in original units if needed
if isfield(model_results, 'Cy')
    model_results.pred_y = (model_results.pred_y .* model_results.Sy) + model_results.Cy;
    model_results.obs_y = (Y .* model_results.Sy) + model_results.Cy;
else
    model_results.obs_y = Y;
end

% Get performance metrics
model_results.r2 = get_model_r2(model_results.obs_y, model_results.pred_y);
model_results.mse = get_mse(model_results.obs_y, model_results.pred_y);

% Save final model if indicated
if cfg.save_model_results == 1
    if isfile('model_results.mat')
       save('model_results.mat', 'mdl_final', '-append');
    else
       save('model_results.mat', 'mdl_final', '-v7.3');
    end
end

end