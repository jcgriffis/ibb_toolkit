function model_results = run_mass_univariate_lr_int(X, Y, cfg, model_results)

% Run mass univariate logistic regression - doesn't currently work
% Joseph Griffis 2024

% Preallocate outputs
pval = zeros(size(X,2),1);
tstat = zeros(size(X,2),1);

% glmfit requires Y to be in range [0,1]
if min(Y) ~= 0 
    Y(Y < 0) = 0;
end

% Get confounds (parfor doesn't like dot indexing)
confounds = cfg.confounds;
int_term = confounds(:,cfg.int_term);

% Define options
options = statset('Display','off', 'MaxIter', 5);

% Parallel processing if indicated
if cfg.parallel == 1
    parfor i = 1:length(tstat)
        warning('off', 'stats:glmfit:IllConditioned');     
        warning('off', 'stats:glmfit:IterationLimit');
        [~,~,stats] = glmfit([X(:,i).*int_term, X(:,i), confounds], Y, 'binomial', 'Options', options);
        pval(i) = stats.p(2);
        tstat(i) = stats.t(2);
    end
else
    warning('off', 'stats:glmfit:IllConditioned');  
    warning('off', 'stats:glmfit:IterationLimit');
    for i = 1:length(tstat)
        [~,~,stats] = glmfit([X(:,i).*int_term, X(:,i), confounds], Y, 'binomial', 'Options', options);
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