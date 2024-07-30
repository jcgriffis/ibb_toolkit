function model_results = run_ols_modeling(X, Y, cfg, model_results)

% Run OLS regression - limited support due to need for cfg.model_terms
% Joseph Griffis 2024

% Create empty field if it doesn't exist
if ~isfield(cfg, 'categorical_X')
    cfg.categorical_X = [];
end

% Fit model to training dataset
if ~isempty(cfg.categorical_X)
    if isfield(cfg, 'model_terms')
        mdl = fitlm(X, Y, cfg.model_terms, "CategoricalVars", cfg.categorical_X);
    else
        mdl = fitlm(X, Y, "CategoricalVars", cfg.categorical_X);
    end
else
    if isfield(cfg, 'model_terms')
        mdl = fitlm(X, Y, cfg.model_terms);
    else
        mdl = fitlm(X, Y);
    end
end

% Store final model
model_results.final_mdl = mdl;

end