function model_results = run_bootstrap_plsda(X, Y, cfg, model_results)

% Run bootstrap analysis for PLSDA
% Joseph Griffis 2024

% Set parallel options
if cfg.parallel == 1
    options = statset('UseParallel', true); % parallel processing
else
    options = statset('UseParallel', false); % no parallel processing
end

%%% Get bootstrapped CIs

% Define bootstrap function
opt_k = model_results.opt_k; % Optimized component number
bootfun = @(IV,DV) run_boot_plsda(IV,DV,opt_k); % returns r-squared, betas, and vip

if cfg.parallel == 1
    % Bootstrap distribution
    bootstat = tall(bootstrp(cfg.boot.n_boot,bootfun,X,Y,'Options',options));    

    % Observed statistics
    stat = tall(bootfun(X,Y)');
else
    % Bootstrap distribution    
    bootstat = bootstrp(cfg.boot.n_boot,bootfun,X,Y,'Options',options);    

    % Observed statistics
    stat = bootfun(X,Y)';  
end

% Get bootstrapped p-values/CIs/etc.
model_results = get_boot_statistics(X, Y, stat, bootstat, bootfun, model_results, cfg);

end