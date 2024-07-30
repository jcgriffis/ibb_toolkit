function model_results = get_boot_cis_cl(stat, bootstat, bootfun, X ,Y, model_results, cfg)

if cfg.parallel == 1
    stat = gather(stat);
    bootstat = gather(bootstat);
end

ci = get_bootstrap_cis(cfg.boot.ci_alpha_thresh_model, stat, bootstat, bootfun, {X,Y}, cfg.boot.ci_type, []);
model_results.boot.classrate_ci = ci(:,1:3); % CI for class rates  
model_results.boot.classrate_ci(model_results.boot.classrate_ci > 1)=1;
model_results.boot.classrate_ci(model_results.boot.classrate_ci < 0)=0;

model_results.boot.roc_auc_ci = ci(:,4); % CI for AUC
clear ci;

% Get CI for model coefficients/feature importances
if cfg.boot.ci_coeff_fwe == 1
    alpha_thresh = cfg.boot.ci_alpha_thresh_coeff ./ size(X,2);
else
    alpha_thresh = cfg.boot.ci_alpha_thresh_coeff;
end

ci = get_bootstrap_cis(alpha_thresh, stat, bootstat, bootfun, {X,Y}, cfg.boot.ci_type, []);  
% threshold betas to retain non-zero crossing CIs
disp('Thresholding betas based on CI widths....');
model_results.boot.coeff_ci = ci(:,5:end); % CIs for betas
model_results.boot.coeff_ci_thresh = get_ci_thresholded_betas(model_results.coeff, ci(:,5:end)); % Thresholded betas

end