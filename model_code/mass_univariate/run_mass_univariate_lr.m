function model_results = run_mass_univariate_lr(X, Y, cfg, model_results)

% Run mass univariate logistic regression - doesn't currently work
% Joseph Griffis 2024

% Preallocate outputs
pval = zeros(size(X,2),1);
tstat = zeros(size(X,2),1);

% glmfit requires Y to be in range [0,1]
if min(Y) == -1 
    Y(Y < 0) = 0;
end

% Parallel processing if indicated
if cfg.parallel == 1
    parfor i = 1:length(tstat)
        [~,~,stats] = glmfit(X(:,i), Y,'binomial', 'Options', statset('Display', 'off'));
        pval(i) = stats.p(2);
        tstat(i) = stats.t(2);
    end
else
    for i = 1:length(tstat)
        [~,~,stats] = glmfit(X(:,i), Y,'binomial', 'Options', statset('Display', 'off'));
        pval(i) = stats.p(2);
        tstat(i) = stats.t(2);
    end
end

% Correct voxel-level p-values 
model_results.coeff_vfwe_pvals = bonf_holm(pval);
model_results.coeff_vfdr_pvals = mafdr(pval, 'BHFDR', true);

% Store results
model_results.coeff_pvals = pval;
model_results.tstat = tstat;
model_results.obs_y = Y;

end