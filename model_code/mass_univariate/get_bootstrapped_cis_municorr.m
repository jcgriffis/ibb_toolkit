function model_results = get_bootstrapped_cis_municorr(X, Y, cfg, model_results)

% Bootstrap CIs for BM test
% Joseph Griffis 2024

% Bootstrap function
bootfun = @(X,Y) corr(X,Y);

% Set parallel options
if cfg.parallel == 1
    options = statset('UseParallel', true);
else
    options = statset('UseParallel', false); % no parallel processing
end

% Run bootstrapping to get bootstat
bootstat = bootstrp(cfg.n_boot,bootfun,X,Y,'Options',options);    

% Get CIs
alpha_thresh = set_bootstrap_alpha_coeff(cfg, size(X,2)); % adjust for multiple predictors if indicated
ci = get_bootstrap_cis(alpha_thresh, bootfun(X,Y)', bootstat, bootfun, {X,Y}, cfg.ci_type, []);
model_results.corr_ci = ci; % CI for R-squared   

% Threshold correlations to remove zero-crossing CIs
cross_zero = find(sign(ci(1,:))-sign(ci(2,:)));
model_results.corr_thresh = model_results.corr;
model_results.corr_thresh(cross_zero)=0;

if cfg.save_boot_results == 1
    model_results.bootstat = bootstat;
end

end

