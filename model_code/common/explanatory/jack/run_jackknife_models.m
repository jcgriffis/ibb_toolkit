function model_results = run_jackknife_models(X, Y, cfg, model_results)

% Run jackknife analysis to get t-statistics and p-values
% Joseph Griffis 2024

% Set parallel options
if cfg.parallel == 1
    options = statset('UseParallel', true); % parallel processing
else
    options = statset('UseParallel', false); % no parallel processing
end

%%% Specify jack-knife analysis parameters based on modeling approach
switch cfg.model_spec
    case 'plsr'
        opt_k = model_results.opt_k; % Optimized component number
        jackfun = @(IV,DV) run_boot_plsr(IV,DV,opt_k); 
    case 'pls_da'
        opt_k = model_results.opt_k; % Optimized component number
        jackfun = @(IV,DV) run_boot_plsda(IV,DV,opt_k); 
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
        jackfun = @(IV, DV) run_boot_regmdl(IV, DV, lambda, learner, reg_type); 
    case {'logistic_ridge', 'logistic_lasso', 'rlinsvc'}
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
        jackfun = @(IV, DV) run_boot_classmdl(IV, DV, lambda, learner, reg_type);
    case {'linsvr', 'kernsvr'}
        C = model_results.C;
        gamma = model_results.gamma;
        epsilon = model_results.epsilon;
        standardize = 0;
        if strcmp(cfg.kernel, 'linear') == 1
            kernel = 0;
        elseif strcmp(cfg.kernel, 'rbf') == 1
            kernel = 1;
        end        
        jackfun = @(IV,DV) run_boot_svr(IV,DV,C,gamma,epsilon,kernel,standardize); 
    case {'linsvc', 'kernsvc'}
        C = model_results.C;
        gamma = model_results.gamma;
        if strcmp(cfg.kernel, 'linear') == 1
            kernel = 0;
        elseif strcmp(cfg.kernel, 'rbf') == 1
            kernel = 1;
        end
        jackfun = @(IV,DV) run_boot_svc(IV,DV,C,gamma,kernel);
end
                
% Get jack-knife estimates
jackstat = jackknife(jackfun,X,Y,'Options',options);

% Drop non-coefficient statistics
if cfg.cat_Y == 0
    jackstat = jackstat(:,2:end); % drop the statistic for r-squared
elseif cfg.cat_Y == 1
    jackstat = jackstat(:,4:end); % drop the statistics for classification accuracy
end

% Compute jack-knife t-statistics and p-values
N = length(Y);
jack_mean = mean(jackstat,1);
model_results.jack.t_coeff = (jack_mean ./ sqrt(((N-1)./N).*sum((jackstat-jack_mean).^2)))';
model_results.jack.t_coeff(isnan(model_results.jack.t_coeff))=0;
model_results.jack.p_coeff = 2 * tcdf(-abs(model_results.jack.t_coeff), N-1);

% Apply multiple comparisons corrections if indicated
if cfg.jack.coeff_fdr == 1
    model_results.jack.fdrp_coeff = mafdr(model_results.jack.p_coeff, 'BHFDR', true);
end
if cfg.jack.coeff_fwe == 1
    model_results.jack.fwep_coeff = bonf_holm(model_results.jack.p_coeff);
end

end