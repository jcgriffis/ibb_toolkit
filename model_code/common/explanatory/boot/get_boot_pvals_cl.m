function model_results = get_boot_pvals_cl(stat, bootstat, model_results, cfg)

% Get bias-corrected version
bias = mean(bootstat-stat,1);

% Get z-values for bootstrapped statistics
if cfg.parallel == 1
    z_stats = gather((stat - bias)./std(bootstat,0,1));
else
    z_stats = (stat - bias)./std(bootstat,0,1);
end

% Convert z-values to p-values
p_boot = 2.*(1-normcdf(abs(z_stats)));

% Separate p-values for different statistics
model_results.boot.p_coeff = p_boot(5:end)';
model_results.boot.z_coeff = z_stats(5:end)';
clear p_boot z_stats

% Do multiple comparisons correction
if cfg.boot.coeff_fdr == 1
    model_results.boot.fdrp_coeff = mafdr(model_results.boot.p_coeff, 'BHFDR', true);
end
if cfg.boot.coeff_fwe == 1
    model_results.boot.fwep_coeff = bonf_holm(model_results.boot.p_coeff);
end

end