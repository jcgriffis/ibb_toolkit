function model_results = run_bootstrap_svr(X, Y, cfg, model_results)

% Run bootstrap analyses
% Joseph Griffis 2024

% Dummy code parameters for bootstrap routine (it doesn't like strings or dot indexing)
C = model_results.C;
gamma = model_results.gamma;
standardize = cfg.standardize;
if strcmp(cfg.kernel, 'linear') == 1
    kernel = 0;
elseif strcmp(cfg.kernel, 'rbf') == 1
    kernel = 1;
end

% define bootstrap function
bootfun = @(IV,DV) run_boot_svr(IV,DV,C,gamma,kernel,standardize); 

if cfg.parallel == 1
    options = statset('UseParallel', true); % parallel processing
else
    options = statset('UseParallel', false); % no parallel processing
end

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
