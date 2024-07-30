function model_results = get_boot_cis_reg(stat, bootstat, bootfun, X,Y, model_results, cfg)

if cfg.parallel == 1
    stat = gather(stat);
    bootstat = gather(bootstat);
end

% Get CI for model R-squared
ci = gather(get_bootstrap_cis(cfg.boot.ci_alpha_thresh_model, stat, bootstat, bootfun, {X,Y}, cfg.boot.ci_type, []));
model_results.boot.r2_ci = ci(:,1); % CI for R-squared    
clear ci;

% Get CI for model coefficients/feature importances
if cfg.boot.ci_coeff_fwe == 1
    alpha_thresh = cfg.boot.ci_alpha_thresh_coeff ./ size(X,2);
else
    alpha_thresh = cfg.boot.ci_alpha_thresh_coeff;
end
ci = gather(get_bootstrap_cis(alpha_thresh, stat, bootstat, bootfun, {X,Y}, cfg.boot.ci_type, []));  

% threshold betas to retain non-zero crossing CIs
disp('Thresholding betas based on CI widths....');
model_results.boot.coeff_ci = ci(:,2:end); % CIs for betas
model_results.boot.coeff_ci_thresh = get_ci_thresholded_betas(model_results.coeff, ci(:,2:end)); % Thresholded betas

end