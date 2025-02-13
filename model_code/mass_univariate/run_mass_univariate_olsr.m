function model_results = run_mass_univariate_olsr(X, Y, cfg, model_results)

% Run mass univariate OLS regression
% Joseph Griffis 2024

% Preallocate outputs
pval = zeros(size(X,2),1);
tstat = zeros(size(X,2),1);

% Get confounds
confounds = cfg.confounds;

% Parallel processing if indicated
if cfg.parallel == 1
    parfor i = 1:length(tstat)
        warning('off', 'stats:glmfit:IllConditioned');     
        warning('off', 'stats:glmfit:IterationLimit');        
        [~,~,stats] = glmfit([X(:,i), confounds], Y,'normal', 'Options', statset('Display', 'off'));
        pval(i) = stats.p(2);
        tstat(i) = stats.t(2);
    end
else
    warning('off', 'stats:glmfit:IllConditioned');     
    warning('off', 'stats:glmfit:IterationLimit');    
    for i = 1:length(tstat)
        [~,~,stats] = glmfit([X(:,i), confounds], Y,'normal', 'Options', statset('Display', 'off'));
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