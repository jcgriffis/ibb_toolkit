function pred_y = predict_new(model_results, X_new)

% Get model predictions for new data
% Joseph Griffis 2024

% Subset out included predictor indices from the new predictor matrix
if size(X_new,2) ~= length(model_results.X_ind)
    X_new = X_new(:, model_results.X_ind);
end

% Apply DTLVC if done for original model
if model_results.cfg.dtlvc == 1
    X_new = apply_dtlvc(X_new);
end

% Use predict function if not a PLS model
if ~strcmp(model_results.cfg.model_spec, 'plsr') && ~strcmp(model_results.cfg.model_spec, 'pls_da')
    
    % Load final model 
    load(fullfile(model_results.cfg.out_dir, 'model_results.mat'), 'mdl_final');

    % Center and scale predictors if needed
    if isfield(model_results, 'Cx')
        X_new = (X_new .* model_results.Sx) + model_results.Cx;
    end        
    
    % Generate new predictions
    if model_results.cfg.cat_Y == 0
        pred_y = predict(mdl_final, X_new);
    else
        [pred_y.labels, pred_y.scores] = predict(mdl_final, X_new);
        pred_y.scores = pred_y.scores(:,2);
        pred_y.labels(:) = -1;
        pred_y.labels(pred_y.scores >= model_results.opt_score_thresh) = 1;
    end

% Otherwise, format the model and new data and apply the model to generate predictions
else
    
    % Center and scale predictors if needed
    if isfield(model_results, 'Cx')
        X_new = (X_new .* model_results.Sx) + model_results.Cx;
    end        
    
    % Generate new predictions
    if model_results.cfg.cat_Y == 0
        pred_y = [ones(size(X_new,1), 1) X_new]*[model_results.beta_0, model_results.coeff']';
    else
        pred_y.scores = [ones(size(X_new,1), 1) X_new]*[model_results.beta_0, model_results.coeff']';        
        pred_y.labels = pred_y.scores;
        pred_y.labels(:) = -1;
        pred_y.labels(pred_y.scores >= model_results.opt_score_thresh) = 1;
    end
end

% Put continuous outcomes back on original scale if original model had standardized Y
if model_results.cfg.cat_Y == 0 && isfield(model_results, 'Cy')
    pred_y = (pred_y .* model_results.Sy) + model_results.Cy;
end

end