function model_results = run_mass_univariate_mnr(X, Y, cfg, model_results)

% Run mass univariate ordinal logistic regression
% Joseph Griffis 2024

% Check version
my_version = version('-release');
my_version = str2double(my_version(1:end-1));

% Preallocate outputs
pval = zeros(size(X,2),1);
tstat = zeros(size(X,2),1);
coeff = zeros(size(X,2),1);

% Get confounds (parfor doesn't like dot indexing)
confounds = cfg.confounds;
n_levels = numel(unique(Y));

% Parallel processing if indicated
if cfg.parallel == 1
    if my_version >= 2023
        parfor i = 1:length(tstat)
            warning('off', 'MATLAB:nearlySingularMatrix');        
            warning('off', 'stats:mnrfit:IterOrEvalLimit');      
            mdl = fitmnr([X(:,i), confounds], Y, 'model', 'ordinal');
            pval(i) = mdl.Coefficients.pValue(n_levels);
            tstat(i) = mdl.Coefficients.tStat(n_levels);
            coeff(i) = mdl.Coefficients.Value(n_levels);
        end
    else
        parfor i = 1:length(tstat)
            warning('off', 'MATLAB:nearlySingularMatrix');        
            warning('off', 'stats:mnrfit:IterOrEvalLimit');      
            [B,~,stats] = mnrfit([X(:,i), confounds], Y, 'model', 'ordinal');
            pval(i) = stats.p(n_levels);
            tstat(i) = stats.t(n_levels);
            coeff(i) = B(n_levels);
        end
    end
else
    if my_version >= 2023
        for i = 1:length(tstst)
            warning('off', 'MATLAB:nearlySingularMatrix');        
            warning('off', 'stats:mnrfit:IterOrEvalLimit');      
            mdl = fitmnr([X(:,i), confounds], Y, 'model', 'ordinal');
            pval(i) = mdl.Coefficients.pValue(n_levels);
            tstat(i) = mdl.Coefficients.tStat(n_levels);
            coeff(i) = mdl.Coefficients.Value(n_levels);
        end
    else
        for i = 1:length(tstat)
            warning('off', 'MATLAB:nearlySingularMatrix');        
            warning('off', 'stats:mnrfit:IterOrEvalLimit');      
            [B,~,stats] = mnrfit([X(:,i), confounds], Y, 'model', 'ordinal');
            pval(i) = stats.p(n_levels);
            tstat(i) = stats.t(n_levels);
            coeff(i) = B(n_levels);
        end
    end
end

% Correct voxel-level p-values 
model_results.coeff_vfwe_pvals = bonf_holm(pval);
model_results.coeff_vfdr_pvals = mafdr(pval, 'BHFDR', true);

% Store results
model_results.coeff_pvals = pval;
model_results.tstat = tstat;
model_results.coeff = coeff;
model_results.obs_y = Y;

end
