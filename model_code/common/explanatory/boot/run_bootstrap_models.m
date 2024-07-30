function model_results = run_bootstrap_models(X, Y, cfg, model_results)

% Run jackknife analysis to get t-statistics and p-values
% Joseph Griffis 2024

% Set parallel options
if cfg.parallel == 1
    options = statset('UseParallel', true); % parallel processing
else
    options = statset('UseParallel', false); % no parallel processing
end

% Return variables to original scale if relevant
if isfield(model_results, 'Cx')
    X = (X .* model_results.Sx) + model_results.Cx;
end      
if isfield(model_results, 'Cy')
    Y = (Y .* model_results.Sy) + model_results.Cy;
end

% Standardization
standardize = cfg.standardize;

%%% Specify bootstrap analysis parameters based on modeling approach
switch cfg.model_spec
    case 'plsr'
        opt_k = model_results.opt_k; % Optimized component number
        bootfun = @(IV,DV) run_boot_plsr(IV,DV,opt_k,standardize); 
    case 'pls_da'
        opt_k = model_results.opt_k; % Optimized component number
        cost = cfg.cost;
        bootfun = @(IV,DV) run_boot_plsda(IV,DV,opt_k,cost,standardize); 
    case {'ridge', 'lasso', 'rlinsvr'}
        reg_type = cfg.reg_type;
        lambda = model_results.Lambda;
        learner = cfg.learner;       
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
        bootfun = @(IV, DV) run_boot_regmdl(IV, DV, lambda, learner, reg_type, standardize); 
    case {'logistic_ridge', 'logistic_lasso', 'rlinsvc'}
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
        bootfun = @(IV, DV) run_boot_classmdl(IV, DV, lambda, learner, reg_type, cost, standardize);
    case {'linsvr', 'kernsvr'}
        C = model_results.C;
        gamma = model_results.gamma;
        if strcmp(cfg.kernel, 'linear') == 1
            kernel = 0;
        elseif strcmp(cfg.kernel, 'rbf') == 1
            kernel = 1;
        end        
        bootfun = @(IV,DV) run_boot_svr(IV,DV,C,gamma,kernel,standardize); 
    case {'linsvc', 'kernsvc'}
        C = model_results.C;
        gamma = model_results.gamma;
        standardize = cfg.standardize;
        cost = cfg.cost;
        if strcmp(cfg.kernel, 'linear') == 1
            kernel = 0;
        elseif strcmp(cfg.kernel, 'rbf') == 1
            kernel = 1;
        end
        bootfun = @(IV,DV) run_boot_svc(IV,DV,C,gamma,kernel,standardize, cost);
end

% Run bootstrap analyses
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