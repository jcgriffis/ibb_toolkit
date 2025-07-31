function [model_results] = run_mass_univariate_corr(X, Y, cfg, model_results)

% Run mass-univariate correlation analyses
% Joseph Griffis 2024

% Run mass univariate correlations
if ~isfield(cfg, 'cor_type')
    cor_type = 'Pearson';
else
    cor_type = cfg.cor_type;
end
if isempty(cfg.confounds)
    disp('No confounds specified...running bivariate correlations')
    [model_results.coeff, model_results.coeff_pvals] = corr(X, Y, 'type', cor_type);
else
    disp('Confounds specified...running partial correlations')
    [model_results.coeff, model_results.coeff_pvals] = partialcorr(X, Y, cfg.confounds, 'type', cor_type);
end

% Correct voxel-level p-values if indicated
model_results.coeff_vfwe_pvals = bonf_holm(model_results.coeff_pvals);
model_results.coeff_vfdr_pvals = mafdr(model_results.coeff_pvals, 'BHFDR', true);

% Get t-statistics
if isempty(cfg.confounds)
    model_results.tstat = model_results.coeff ./ sqrt((1-model_results.coeff.^2) ./ (length(Y)-2));
else
    model_results.tstat = model_results.coeff .* sqrt((length(Y)-2-size(cfg.confounds,2)) ./ (1-model_results.coeff.^2));
end

if isempty(cfg.confounds)
    % Generate CPM model
    model_results.obs_y = Y;
    cpm_train = X*model_results.coeff;
    if numel(unique(Y)) > 2
        mdl = fitlm(cpm_train, Y);
    else
        mdl = fitlm(Y, cpm_train, 'CategoricalVars', 'x1');
    end
    model_results.cpm_mdl = mdl;
end

end