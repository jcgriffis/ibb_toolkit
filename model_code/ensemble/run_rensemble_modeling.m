function [model_results] = run_rensemble_modeling(X, Y, cfg, model_results)

% Use K-fold cross-validation to tune classification ensemble hyper-parameters and fit model
% Joseph Griffis 2024
% (heavily copied from run_censemble_modeling.m by Stein Acker)

if cfg.optimize_hyperparams == 1

    % Get hyper-parameter options and cross-validation partitions
    model_results.hp_opt = get_hp_opts(cfg, length(Y));

    % Run optimization
    if iscellstr(cfg.hp_opt.to_optimize)
        if contains('Method', cfg.hp_opt.to_optimize)
            mdl_final = fitrensemble(X,Y,'Learners', 'tree',...
                'OptimizeHyperparameters', cfg.hp_opt.to_optimize,...
                'HyperparameterOptimizationOptions', model_results.hp_opt);

        else
            mdl_final = fitrensemble(X,Y,'Learners', 'tree',...
                'Method', cfg.method,...
                'Options', statset('UseParallel',cfg.parallel));
        end
    elseif isa(cfg.hp_opt.to_optimize(1), 'optimizableVariable')
        if cfg.hp_opt.to_optimize(1).Optimize == 1
            mdl_final = fitrensemble(X,Y,'Learners', 'tree',...
                'OptimizeHyperparameters', cfg.hp_opt.to_optimize,...
                'HyperparameterOptimizationOptions', model_results.hp_opt);
        else
            mdl_final = fitrensemble(X,Y,'Learners', 'tree',...
                'Method',cfg.method,...
                'OptimizeHyperparameters', cfg.hp_opt.to_optimize,...
                'HyperparameterOptimizationOptions', model_results.hp_opt);
        end
    end
else
    mdl_final = fitrensemble(X,Y,'Learners', 'tree',...
        'Method', cfg.method,...
        'Options', statset('UseParallel',cfg.parallel));
end

% Get optimal hyper-parameters
model_results.method = mdl_final.ModelParameters.Method;
model_results.n_learn = mdl_final.ModelParameters.NLearn;

% Get fitted Y
model_results.pred_y = predict(mdl_final, X, 'ObservationsIn', 'rows');

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