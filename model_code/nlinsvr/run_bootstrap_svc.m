function model_results = run_bootstrap_svc(X, Y, cfg, model_results)

% Run bootstrap analyses
% Joseph Griffis 2024

% Dummy code parameters for bootstrap routine (it doesn't like strings or dot indexing)
C = model_results.C;
gamma = model_results.gamma;
standardize = cfg.standardize;
cost = cfg.cost;
if strcmp(cfg.kernel, 'linear') == 1
    kernel = 0;
elseif strcmp(cfg.kernel, 'rbf') == 1
    kernel = 1;
end

% define bootstrap function
bootfun = @(IV,DV) run_boot_svc(IV,DV,C,gamma,kernel,standardize,cost); % returns R^2, AIC, and betas as a vector

if cfg.parallel == 1
    options = statset('UseParallel', true); % parallel processing
else
    options = statset('UseParallel', false); % no parallel processing
end

% Get bootstrap weights
w = zeros(size(Y));
p1 = numel(Y(Y==1))./numel(Y); % prevalence of class 1
p2 = numel(Y(Y==-1))./numel(Y); % prevalence of class 2 
w1 = (1./length(Y)) ./ p1; % weights for class 1 
w2 = (1./length(Y)) ./ p2; % weights for class 2 
w(Y==1) = w1./2; % divide by 2 so weights sum to 1 
w(Y==-1)= w2./2;

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