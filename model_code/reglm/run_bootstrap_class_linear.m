function model_results = run_bootstrap_class_linear(X, Y, cfg, model_results)

% Run bootstrap analyses
% Joseph Griffis 2024

% Dummy code learner for bootstrap routine (it doesn't like strings)
reg_type = cfg.reg_type;
lambda = model_results.Lambda;
learner = cfg.learner;
cost = cfg.cost;

if strcmp(learner, 'svm') == 1
    learner = 1;
else
    learner = 0;
end
if strcmp(reg_type, 'lasso') == 1
    reg_type = 1;
else
    reg_type = 0;
end
% define bootstrap function
bootfun = @(IV, DV) run_boot_classmdl(IV, DV, lambda, learner, reg_type, cost); % returns R^2, AIC, and betas as a vector
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
