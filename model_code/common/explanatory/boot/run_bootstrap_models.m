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

%%% Specify bootstrap analysis parameters based on modeling approach
switch cfg.model_spec
    case 'plsr'
        bootfun = @(IV,DV) run_boot_plsr(IV,DV,model_results.opt_k,cfg.standardize, cfg.standardize_method, cfg.standardize_type); 
    case 'pls_da'
        bootfun = @(IV,DV) run_boot_plsda(IV,DV, model_results.opt_k, cfg.cost, cfg.standardize, cfg.standardize_method, cfg.standardize_type); 
    case {'ridge', 'lasso', 'rlinsvr'}
        bootfun = @(IV, DV) run_boot_regmdl(IV, DV, model_results.Lambda, cfg.learner, cfg.reg_type, cfg.standardize, cfg.standardize_method, cfg.standardize_type); 
    case {'logistic_ridge', 'logistic_lasso', 'rlinsvc'}
        bootfun = @(IV, DV) run_boot_classmdl(IV, DV, model_results.lambda, cfg.learner, cfg.reg_type, cfg.cost, cfg.standardize, cfg.standardize_method, cfg.standardize_type);
    case {'linsvr', 'kernsvr'}
        bootfun = @(IV,DV) run_boot_svr(IV,DV,model_results.C,model_results.gamma,cfg.kernel,model_results.epsilon,cfg.standardize,cfg.standardize_method, cfg.standardize_type); 
    case {'linsvc', 'kernsvc'}
        bootfun = @(IV,DV) run_boot_svc(IV,DV,model_results.C,model_results.gamma,cfg.kernel,cfg.standardize,cfg.cost,cfg.standardize_method,cfg.standardize_type);
end

% Weight observations by stratification groups if applicable
if isfield(cfg, 'strat_groups') && isfield(cfg.boot, 'weighted')

    if cfg.boot.weighted == 1

        % Get bootstrap weights
        strat_groups = unique(cfg.strat_groups);
        w = zeros(size(cfg.strat_groups));
    
        for i = 1:length(strat_groups)
    
            p = numel(cfg.strat_groups(cfg.strat_groups==strat_groups(i)))./numel(cfg.strat_groups); % prevalence of group
            w1 = (1./length(cfg.strat_groups)) ./ p; % weights for group
            w(cfg.strat_groups==strat_groups(i)) = w1./length(strat_groups); % divide by N so weights sum to 1 
    
        end
        
    else

        % Default weighting
        w = ones(length(Y),1)./length(Y);        
        
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