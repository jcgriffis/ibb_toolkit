function model_results = run_mass_univariate_mnr(X, Y, cfg, model_results)

% Run mass univariate logistic regression - doesn't currently work
% Joseph Griffis 2024

% Preallocate outputs
pval = zeros(size(X,2),1);
tstat = zeros(size(X,2),1);

% Get confounds (parfor doesn't like dot indexing)
confounds = cfg.confounds;
n_levels = numel(unique(Y));

% Parallel processing if indicated
if cfg.parallel == 1
    parfor i = 1:length(tstat)
        warning('off', 'MATLAB:nearlySingularMatrix');        
        warning('off', 'stats:mnrfit:IterOrEvalLimit');      
        [~,~,stats] = mnrfit([X(:,i), confounds], Y, 'model', 'ordinal');
        pval(i) = stats.p(n_levels);
        tstat(i) = stats.t(n_levels);
    end
else
    for i = 1:length(tstat)
        warning('off', 'MATLAB:nearlySingularMatrix');        
        warning('off', 'stats:mnrfit:IterOrEvalLimit');
        [~,~,stats] = mnrfit([X(:,i), confounds], Y, 'model', 'ordinal');
        pval(i) = stats.p(n_levels);
        tstat(i) = stats.t(n_levels);
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