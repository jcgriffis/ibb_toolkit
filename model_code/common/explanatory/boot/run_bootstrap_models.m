function model_results = run_bootstrap_models(X, Y, cfg, model_results)

% Run bootstrap analysis to get t-statistics and p-values
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
if cfg.standardize > 0 && ~isfield(cfg, 'standardize_method')
    cfg.standardize_method = 'zscore';
    cfg.standardize_type = 'std';
elseif cfg.standardize > 0 && isfield(cfg, 'standardize_method')
    if strcmp(cfg.standardize_method, 'range') && ~isfield(cfg, 'standardize_type')
        cfg.standardize_type = [0,1];
    elseif strcmp(cfg.standardize_method, 'zscore') && ~isfield(cfg, 'standardize_type')
        cfg.standardize_type = 'std';
    end
elseif cfg.standardize == 0 
    cfg.standardize_method = [];
    cfg.standardize_type = [];
end
method = cfg.standardize_method;
type = cfg.standardize_type;

%%% Specify bootstrap analysis parameters based on modeling approach
switch cfg.model_spec
    case 'plsr'
        opt_k = model_results.opt_k; % Optimized component number
        bootfun = @(IV,DV) run_boot_plsr(IV,DV,opt_k,standardize, method, type); 
    case 'pls_da'
        opt_k = model_results.opt_k; % Optimized component number
        cost = cfg.cost;
        bootfun = @(IV,DV) run_boot_plsda(IV,DV,opt_k,cost,standardize, method, type); 
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
        bootfun = @(IV, DV) run_boot_regmdl(IV, DV, lambda, learner, reg_type, standardize, method, type); 
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
        bootfun = @(IV, DV) run_boot_classmdl(IV, DV, lambda, learner, reg_type, cost, standardize, method, type);
    case {'linsvr', 'kernsvr'}
        C = model_results.C;
        gamma = model_results.gamma;
        if strcmp(cfg.kernel, 'linear') == 1
            kernel = 0;
        elseif strcmp(cfg.kernel, 'rbf') == 1
            kernel = 1;
        end        
        bootfun = @(IV,DV) run_boot_svr(IV,DV,C,gamma,kernel,standardize, method, type); 
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
        bootfun = @(IV,DV) run_boot_svc(IV,DV,C,gamma,kernel,standardize, cost, method, type);
end

% Weight observations by stratification groups if applicable
if isfield(cfg, 'strat_groups')

    % Get bootstrap weights
    strat_groups = unique(cfg.strat_groups);
    w = zeros(size(cfg.strat_groups));

    for i = 1:length(strat_groups)

        p = numel(cfg.strat_groups(cfg.strat_groups==strat_groups(i)))./numel(cfg.strat_groups); % prevalence of group
        w1 = (1./length(cfg.strat_groups)) ./ p; % weights for group
        w(cfg.strat_groups==strat_groups(i)) = w1./length(strat_groups); % divide by 2 so weights sum to 1 

    end

else

    % Default weighting
    w = ones(length(Y),1)./length(Y);

end

% Run bootstrap analyses
if cfg.parallel == 1
    % Bootstrap distribution
    bootstat = tall(bootstrp(cfg.boot.n_boot,bootfun,X,Y,'Options',options, 'weights', w));    

    % Observed statistics
    stat = tall(bootfun(X,Y)');
else
    % Bootstrap distribution    
    bootstat = bootstrp(cfg.boot.n_boot,bootfun,X,Y,'Options',options, 'weights', w);    

    % Observed statistics
    stat = bootfun(X,Y)';  
end

% Get bootstrapped p-values/CIs/etc.
model_results = get_boot_statistics(X, Y, stat, bootstat, bootfun, model_results, cfg);

end