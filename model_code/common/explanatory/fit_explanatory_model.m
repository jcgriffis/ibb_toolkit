function model_results = fit_explanatory_model(X, Y, model_results, cfg)

% Fit explanatory model to full dataset
% Performs confound regression if specified in cfg
% Standardizes data if specified in cfg
% Performs permutation testing if specified in cfg
% Performs bootstrap analyses if specified in cfg

% Joseph Griffis 2024

%%%%% Regress out confound variables if they're included
if ~isempty(cfg.confounds)
    disp('Regressing confounds out of Y...');    
    [Y, model_results.r2_confound, model_results.coeff_confound] = regress_confounds_from_Y(Y, cfg.confounds);
    disp(['Confound model explained ' num2str(round(model_results.r2_confound, 4)*100) '% variance in Y.']);
    disp('Continuing analysis with Y residuals as target variable');
end

% Standardize data as indicated in cfg
if cfg.standardize == 1 % Convert X to z-scores
    if cfg.cat_Y == 0
        [X, model_results.Cx, model_results.Sx] = normalize(X);
        X(isnan(X))=0; % Since normalize will cause columns to become NaN if they are all 0
    elseif cfg.cat_Y == 1
        [X, model_results.Cx, model_results.Sx] = normalize(X);
        X(isnan(X))=0; % Since normalize will cause columns to become NaN if they are all 0
    end
elseif cfg.standardize == 2 % Convert Y to z-scores 
    if cfg.cat_Y == 0
        [Y, model_results.Cy, model_results.Sy] = normalize(Y);
    else
        disp('Cannot standardize categorical outcome - ignoring standardization flag and seting to 0.');
        cfg.standardize = 0;
    end
elseif cfg.standardize == 3 % Convert X and Y to z-scores
    [X, model_results.Cx, model_results.Sx] = normalize(X);
    X(isnan(X))=0; % Since normalize will cause columns to become NaN if they are all 0
    if cfg.cat_Y == 0
        [Y, model_results.Cy, model_results.Sy] = normalize(Y);    
    else
        cfg.standardize = 1;
        disp('Cannot standardize categorical outcome - ignoring standardization flag and seting to 1.');
    end        
end

% Fit full-sample model if indicated
switch cfg.model_spec 
    case 'plsr'
        model_results = run_plsr_modeling(X, Y, cfg, model_results);
    case 'pls_da'
        model_results = run_plsda_modeling(X, Y, cfg, model_results);   
    case {'ridge', 'lasso', 'rlinsvr'} 
        model_results = run_reg_linear_modeling(X, Y, cfg, model_results);
    case {'logistic_ridge', 'logistic_lasso', 'rlinsvc'}
        model_results = run_class_linear_modeling(X, Y, cfg, model_results);
    case {'kernsvr', 'linsvr'}
        model_results = run_svr_modeling(X, Y, cfg, model_results);
    case {'kernsvc', 'linsvc'}
        model_results = run_svc_modeling(X, Y, cfg, model_results);
    case {'censemble'}
        model_results = run_censemble_modeling(X, Y, cfg, model_results);
    case {'rensemble'}
        model_results = run_rensemble_modeling(X, Y, cfg, model_results);
    case 'olsr'
        model_results = run_ols_modeling(X, Y, cfg, model_results);
    case 'municorr'
        model_results = run_mass_univariate_corr(X, Y, model_results);
    case 'munilr'
        model_results = run_mass_univariate_lr(X, Y, cfg, model_results);
    case 'bmunz'
        model_results = run_brunner_munzel_test(X, Y, cfg, model_results);
    case 'ttest'
        model_results = run_ttests(X, Y, cfg, model_results);
end

% Permutation tests
cd(cfg.out_dir);
if cfg.permutation == 1
    switch cfg.model_spec
        case 'plsr'
            [model_results] = run_perm_plsr(X, Y, cfg, model_results);
        case 'pls_da'
            [model_results] = run_perm_plsda(X, Y, cfg, model_results);
        case {'ridge', 'lasso', 'rlinsvr'} 
            [model_results] = run_perm_reglm(X, Y, cfg, model_results);
        case {'logistic_ridge', 'logistic_lasso', 'rlinsvc'}
            [model_results] = run_perm_classlm(X, Y, cfg, model_results);
        case {'kernsvr', 'linsvr'}
            [model_results] = run_perm_svr(X, Y, cfg, model_results);
        case {'kernsvc', 'linsvc'}
            [model_results] = run_perm_svc(X, Y, cfg, model_results);  
        case 'municorr'
            [model_results] = run_perm_mass_uni_corr(X, Y, cfg, model_results);
        case 'munilr'
            [model_results] = run_perm_mass_uni_lr(X, Y, cfg, model_results);  
    end    
end

%%% Bootstrap z-statistics and CIs
if cfg.bootstrap == 1
    disp(['Running bootstrap analyses with ' num2str(cfg.boot.n_boot) ' bootstrap resamples...'])
    model_results = run_bootstrap_models(X, Y, cfg, model_results);
end

%%% Jack-knife t-statistics
if cfg.jackknife == 1
    disp('Running jack-knife analyses to get coefficient test statistics and p-values...')
    model_results = run_jackknife_models(X, Y, cfg, model_results);
end

end