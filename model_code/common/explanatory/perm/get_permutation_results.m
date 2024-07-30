function model_results = get_permutation_results(model_results, perm_coeff, cfg)

% Get coefficient p-values and apply corrections if indicated

% Joseph Griffis 2024

% Compute whole-brain continuous FWE-corrected p-values
disp('Computing permutation p-values...')
if cfg.perm.coeff_cfwe == 1
    disp('Computing whole-brain cFWE p-value thresholds...')
    model_results = get_permutation_coefficient_pvalues(model_results, 'cfwe', perm_coeff, cfg.fwe_thresh, cfg.perm.n_perm);
end

% Compute voxel-level p-values
if cfg.perm.coeff_p == 1
    disp('Computing voxel-level permutation p-values...')
    if ~strcmp(cfg.model_spec, 'ttest')
        model_results.perm.p_coeff = get_perm_beta_pvals(model_results.coeff, perm_coeff, cfg.perm.n_perm);            
    else
        model_results.perm.p_coeff = get_perm_beta_pvals(model_results.tstat, perm_coeff, cfg.perm.n_perm);            
    end        
    if cfg.perm.coeff_fwe == 1
        model_results.perm.fwep_coeff = bonf_holm(model_results.perm.p_coeff);
    end
    if cfg.perm.coeff_fdr == 1
        model_results.perm.fdrp_coeff = mafdr(model_results.perm.p_coeff, 'BHFDR', true);
    end
end

end