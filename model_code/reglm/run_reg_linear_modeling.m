function [model_results] = run_reg_linear_modeling(X, Y, cfg, model_results)

% Use K-fold cross-validation to tune regularized regression hyper-parameters
% and fit final model to full dataset. 
% Joseph Griffis 2023

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
                
                disp(['Hyper-parameter optimization iteration:' num2str(i) '/' num2str(repeats)]) 
                
                mdl = fitrlinear(X,Y,...
                'ObservationsIn', 'rows',...
                'Regularization', reg_type,...
                'Learner', learner,...
                'Lambda', lambda_range,...
                'CVPartition', my_part,...
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

                mdl = fitrlinear(X,Y,...
                'ObservationsIn', 'rows',...
                'Regularization', reg_type,...
                'Learner', learner,...
                'Lambda', lambda_range,...
                'CVPartition', cv_part,...
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
        mdl_final = fitrlinear(X,Y,...
        'ObservationsIn', 'rows',...
        'Regularization', cfg.reg_type,...
        'Lambda', opt_l,...
        'Learner', cfg.learner,...
        'BetaTolerance', 0);
           
    else
            disp('Optimizing lambda with Bayesian optmization');
            
            model_results.hp_opt = get_regmdl_hp_opts(cfg, length(Y));
      
            mdl_final = fitrlinear(X,Y,...
            'ObservationsIn', 'rows',...
            'Regularization', cfg.reg_type,...
            'OptimizeHyperparameters', {'Lambda'},...
            'Learner', cfg.learner,...
            'HyperparameterOptimizationOptions', model_results.hp_opt,...
            'BetaTolerance', 0);
    
    end
else
    if ~isfield(cfg, 'lambda')
        mdl_final = fitrlinear(X,Y,...
        'ObservationsIn', 'rows',...
        'Regularization', cfg.reg_type,...
        'Learner', cfg.learner,...
        'BetaTolerance', 0);   
    else
        mdl_final = fitrlinear(X,Y,...
        'ObservationsIn', 'rows',...
        'Regularization', cfg.reg_type,...
        'Learner', cfg.learner,...
        'BetaTolerance', 0,...
        'Lambda', cfg.lambda);
    end
end

% Get final model R-squared
model_results.pred_y = predict(mdl_final, X, 'ObservationsIn', 'rows');

% Put predictions and observations back in original units if needed
if isfield(model_results, 'Cy')
    model_results.pred_y = (model_results.pred_y .* model_results.Sy) + model_results.Cy;
    model_results.obs_y = (Y .* model_results.Sy) + model_results.Cy;
else
    model_results.obs_y = Y;
end

% Get model performance metrics
model_results.r2 = get_model_r2(model_results.obs_y, model_results.pred_y);
model_results.mse = get_mse(model_results.obs_y, model_results.pred_y);

% Get beta weights
model_results.coeff = mdl_final.Beta;
model_results.beta0 = mdl_final.Bias;
model_results.Lambda = opt_l;

% Save final model if indicated
if cfg.save_model_results == 1
    if isfile('model_results.mat')
       save('model_results.mat', 'mdl_final', '-append');
    else
       save('model_results.mat', 'mdl_final', '-v7.3');
    end
end

end