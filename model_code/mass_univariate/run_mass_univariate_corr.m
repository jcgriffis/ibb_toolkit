function [model_results] = run_mass_univariate_corr(X, Y, model_results)

% Run mass-univariate correlation analyses
% Joseph Griffis 2024

% Run mass univariate correlations
[model_results.coeff, model_results.coeff_pvals] = corr(X, Y);

% Correct voxel-level p-values if indicated
model_results.coeff_vfwe_pvals = bonf_holm(model_results.coeff_pvals);
model_results.coeff_vfdr_pvals = mafdr(model_results.coeff_pvals, 'BHFDR', true);

% Get t-statistics
model_results.tstat = model_results.coeff ./ sqrt((1-model_results.coeff.^2) ./ (length(Y)-2));

% Generate CPM model
model_results.obs_y = Y;
cpm_train = X*model_results.coeff;
mdl = fitlm(cpm_train, Y);
model_results.cpm_mdl = mdl;

end