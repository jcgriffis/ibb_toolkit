function model_results = get_bootstrapped_cis_munilr(X, Y, cfg, model_results)

% Run mass univariate logistic regressions - doesn't currently work
% Joseph Griffis 2024

% Set Y to range for glmfit
Y = Y > 0;

% Set parallel options
if cfg.parallel == 1
    options = statset('UseParallel', true);
    parallel = 1;
else
    options = statset('UseParallel', false); % no parallel processing
    parallel = 0;
end

% Bootstrap function
bootfun = @(X,Y) boot_munilr(X,Y,parallel);

% Run bootstrapping to get bootstat
model_results.bootstat = bootstrp(cfg.n_boot,bootfun,X,Y,'Options',options);    

% Get CIs
alpha_thresh = set_bootstrap_alpha_coeff(cfg, size(X,2)); % adjust for multiple predictors if indicated
ci = get_bootstrap_cis(alpha_thresh, bootfun(X,Y)', model_results.bootstat, bootfun, {X,Y}, cfg.ci_type, []);
model_results.beta_ci = ci; 

% Threshold correlations to remove zero-crossing CIs
cross_zero = find(sign(ci(1,:))-sign(ci(2,:)));
model_results.beta_thresh = model_results.betas;
model_results.beta_thresh(cross_zero)=0;

end

function b = boot_munilr(X,Y,parallel)

% Preallocate
b = zeros(size(X,2),1);
if parallel == 1
    parfor i = 1:length(b)
        [betas,~,~] = glmfit(X(:,i), Y,'binomial');
        b(i) = betas(2);
    end
else
    for i = 1:length(b)
        [betas,~,~] = glmfit(X(:,i), Y,'binomial');
        b(i) = betas(2);
    end
end

end