function cfg = get_default_model_cfg(cfg)

% Generates a pre-filled cfg structure containing modeling parameters and analysis options
% Joseph Griffis 2024

% General analysis fields that should be provided to cfg before inputting to this function (values should be 1 or 0):
%
% model_spec - model to run (value should be one of the values output from calling get_model_option_list()).
% fit_explanatory_model - run inferential modeling on full dataset (value should be 1 or 0)
% optimize_hyperparams - optimze hyper-parameters for multivariate models (value should be 1 or 0)
% cross_validation - run cross-validation analyses (value should be 1 or 0)
% bootstrap - run bootstrap analyses on inferential model (value should be 1 or 0)
% permutation - run permutation tests on inferential model (value should be 1 or 0)
% jackknife - run jack-knife tests on inferential model (value should be 1 or 0)

% Add common directory to path
restore_model_path(cfg.model_dir);
addpath(genpath(fullfile(cfg.model_dir, 'common')));

% Set defaults for general analysis fields if not set
cfg = check_and_fill_cfg_fields(cfg);

% Model-specific options
switch cfg.model_spec
    case {'plsr', 'pls_da'} % PLS models
        cfg.hp_opt.comp_method = 'dif_mse'; % identifies last component that minimized MSE before first component that caused MSE increase (standard method)
        cfg.hp_opt.n_comp = 5; % Max number of model components        
        addpath(fullfile(cfg.model_dir, 'plsr'));
        cfg.hp_opt.bayes_opt = 0;
        if strcmp(cfg.model_spec, 'pls_da')
            cfg.cost = [0, 1; 1, 0];
            cfg.cat_Y = 1;
        else
            cfg.cat_Y = 0;
        end        
        cfg.standardize = 0;
    case {'lasso', 'ridge', 'rlinsvr', 'logistic_lasso', 'logistic_ridge', 'rlinsvc'} % Regularized linear models
        if contains(cfg.model_spec, {'lasso', 'ridge'})
            cfg.learner = 'leastsquares'; % OLS regression
        elseif contains(cfg.model_spec, {'logistic'})
            cfg.learner = 'logistic';
        else
            cfg.learner = 'svm'; % Support vector regression
        end
        if contains(cfg.model_spec, {'logistic', 'svc'})
            cfg.cost = [0, 1; 1, 0];
            cfg.cat_Y = 1;
        else
            cfg.cat_Y = 0;
        end
        if contains(cfg.model_spec, 'lasso')
            cfg.reg_type = 'lasso'; % L1 lasso regularization
        else
            cfg.reg_type = 'ridge'; % L2 ridge regularization
        end
        cfg.hp_opt.bayes_opt = 0;     
        addpath(fullfile(cfg.model_dir, 'reglm'));  
        if ~contains(cfg.model_spec, {'svc', 'logistic'})
            cfg.standardize = 3;
        else
            cfg.standardize = 1;
        end
        cfg.standardize_method = 'range';        
        cfg.standardize_type = [0, 1];
    case {'linsvr', 'kernsvr', 'linsvc', 'kernsvc'} % Support vector machines
        if contains(cfg.model_spec, 'lin')
            cfg.kernel = 'linear'; % linear SVR kernel
        else
            cfg.kernel = 'rbf';
        end
        if contains(cfg.model_spec, 'svc')
            cfg.cost = [0, 1; 1, 0];
            cfg.cat_Y = 1;
        else
            cfg.cat_Y = 0;
        end
        % Optimization parameters
        cfg.hp_opt.box_constraint.optimize = 1; % Box constraint
        cfg.hp_opt.box_constraint.range = [0.001, 100];
        cfg.hp_opt.kernel_scale.optimize = 1; % Kernel scale
        cfg.hp_opt.kernel_scale.range = [0.001, 100];
        cfg.hp_opt.kernel_function.optimize = 0; % Kernel function (linear, polynomial, rbf)
        cfg.hp_opt.poly_order.optimize = 0; % Polynomial order
        cfg.hp_opt.standardize.optimize = 0; % Standardization (separate from toolkit standardization)
        cfg.hp_opt.bayes_opt = 1; % Do Bayesian hyper-parameter optimization
        cfg.hp_opt.opt_iter = 100; % Number of objective function evaluations
        cfg.hp_opt.repartition = true; % Repartition into train/test at each iteration         
        addpath(fullfile(cfg.model_dir, 'nlinsvr'));
        if contains(cfg.model_spec, 'svc')
            cfg.standardize = 1;
            cfg.standardize_method = 'range';
            cfg.standardize_type = [0, 1];        
        else
            cfg.standardize = 3;
            cfg.standardize_method = 'range';
            cfg.standardize_type = [0, 1];            
        end
    case {'rensemble', 'censemble'} % Ensemble predictive models
        if contains(cfg.model_spec, 'censemble')
            cfg.cost = [0, 1; 1, 0];
            cfg.cat_Y = 1;
        else
            cfg.cat_Y = 0;
        end
        cfg.learners = "tree";
        if ~isfield(cfg, 'method') || (isfield(cfg, 'method') && strcmp(cfg.method, 'Bag'))
            cfg.method = 'Bag';
            if strcmp(cfg.model_spec, 'rensemble')
                params = hyperparameters('fitrensemble', [], [], 'tree');
            else
                params = hyperparameters('fitcensemble', [], [], 'tree');
            end
            params(2).Range = [10, 200];
            params(1).Optimize = false;
            params(3).Optimize = false; % Can't optimize LearnRate for Bagged ensemble
            cfg.hp_opt.to_optimize = params;
        elseif isfield(cfg, 'method') && contains(cfg.method, {'LSBoost', 'AdaBoostM1', 'AdaBoostM2', 'GentleBoost', 'LogitBoost'})
            if strcmp(cfg.model_spec, 'rensemble')
                params = hyperparameters('fitrensemble', [], [], 'tree');
            else
                params = hyperparameters('fitcensemble', [], [], 'tree');
            end
            params(2).Range = [10, 200];
            params(1).Optimize = false;
            params(3).Optimize = true; 
            cfg.hp_opt.to_optimize = params;
        else
            cfg.hp_opt.to_optimize = 'auto';
        end
        cfg.hp_opt.bayes_opt = 1; % Do Bayesian hyper-parameter optimization
        cfg.hp_opt.opt_iter = 60; % Number of objective function evaluations
        cfg.hp_opt.repartition = true; % Repartition into train/test at each iteration
        cfg.standardize = 0;
        addpath(fullfile(cfg.model_dir, 'ensemble'));         
    case 'olsr' % Ordinary least squares regression - incomplete and implemented for limited use case (continuous outcome)
        cfg.cat_Y = 0;
        cfg.hp_opt = [];
        cfg.standardize = 0;        
        addpath(fullfile(cfg.model_dir, 'mass_univariate'));
    case 'municorr' % Mass univariate correlation - bivariate Pearson correlations (also generates t-statistics)
        cfg.cat_Y = 0;
        cfg.hp_opt = [];        
        cfg.standardize = 0;        
        addpath(fullfile(cfg.model_dir, 'mass_univariate'));      
    case 'munilr' % Mass univariate logistic regression
        cfg.cat_Y = 1;
        cfg.hp_opt = [];        
        cfg.standardize = 0;        
        addpath(fullfile(cfg.model_dir, 'mass_univariate'));    
    case 'munimnr' % Mass univariate ordinal regression 
        cfg.cat_Y = 0;
        cfg.hp_opt = [];
        cfg.standardize = 0;
        addpath(fullfile(cfg.model_dir, 'mass_univariate'));            
    case 'muniolsr'
        cfg.cat_Y = 0;
        cfg.hp_opt = [];
        cfg.standardize = 0;        
        addpath(fullfile(cfg.model_dir, 'mass_univariate'));
    case 'bmunz' % Brunner-Munzel test - limited functionality (only for fit_explanatory_model, no bootstrap functionality)
        cfg.cat_Y = 0;
        cfg.hp_opt = [];        
        cfg.standardize = 0;        
        addpath(fullfile(cfg.model_dir, 'mass_univariate'));
    case 'ttest' % Mass univariate t-tests - limited functionality (only for fit_explanatory_model, no bootstrap functionality)
        cfg.cat_Y = 0;
        cfg.hp_opt = [];        
        cfg.standardize = 0;        
        cfg.vartype = 'unequal'; % unequal variances t-test
        addpath(fullfile(cfg.model_dir, 'mass_univariate'))
    case 'prop_sub' % Proportional subtraction analysis
        cfg.cat_Y = 1;
        cfg.hp_opt = [];
        cfg.standardize = 0;        
        addpath(fullfile(cfg.model_dir, 'mass_univariate'));
    otherwise
        warning('Model type not provided as input - make sure to set model-specific CFG fields since they will not be set by default!' );
end

% Analysis-specific options
switch cfg.modality
    case 'lesion' 
        disp('Analysis type is lesion - setting defaults')
        if ~isfield(cfg, 'mask_path')
            error('Brain mask required for lesion analysis - please supply brain mask in same space/orientation/dimensions as lesion data')
        end
        
        % Load brain mask
        mask = niftiread(cfg.mask_path);
        mask(isnan(mask))=0;

        cfg.mask_inds = find(mask); % get map indices for non-zero voxels
        cfg.lsm = 1; % lesion-symptom mapping analysis
        cfg.trim_X = 1; % trim X to only include predictors with sufficient observations
        cfg.min_obs = 1; % minimum value to consider for an observation
        cfg.freq_thresh = 10; % minimum frequency of observations with values > min_obs for columns to include when trimming
        cfg.gen_freq_map = 1; % generate lesion frequency map
        cfg.gen_lvol_map = 1; % Generate lesion volume correlation map
    otherwise
        disp('Analysis type is other - setting defaults for non-lesion analysis')
        cfg.lsm = 0;
        if isfield(cfg, 'mask_path')
            mask = niftiread(cfg.mask_path);
            mask(isnan(mask))=0;
            cfg.mask_inds = find(mask);
        end
        cfg.trim_X = 1;
        cfg.min_obs = 0.01;
        cfg.freq_thresh = 10;
        cfg.gen_freq_map = 0;    
        cfg.gen_lvol_map = 0;
end

%%%% Misc options
if ~isfield(cfg, 'dtlvc')
    cfg.dtlvc = 0; % direct total lesion volume control (Zhang et al., 2014)
end
if ~isfield(cfg, 'strat_var')
    cfg.strat_var = []; % stratification variable (used to stratify groups for CV)
    cfg.strat_name = []; % name of stratification variable
end
if ~isfield(cfg, 'confounds')
    cfg.confounds = []; % confound variables (will be regressed out of outcome)
end

%%%% Modeling options

% Within-Dataset Explanatory modeling (or fitting data to full dataset to predict on independent held-out dataset)
if cfg.fit_explanatory_model == 1

    % Jackknife options for explanatory model
    if cfg.jackknife == 1
        
        cfg.jack.write_jack_images = 1; % save images of statistics
        cfg.jack.write_uncorrected_images = 1; % save out uncorrected statistical maps
        cfg.jack.get_pvals = 1; % get p-values for model coefficients
        cfg.jack.coeff_fwe = 1; % apply family-wise error correction to coefficient p-values (Bonferroni-Holm Method)
        cfg.jack.coeff_fdr = 1; % apply false discovery rate correction to coefficient p-values (Benjamini-Hochberg Method)

    end            

    % Bootstrap options for explanatory model
    if cfg.bootstrap == 1

        cfg.boot.n_boot = 1000; % number of bootstrap iterations (minimum should be 1000, 5000 preferred)
        cfg.boot.save_boot_results = 0; % save out bootstrap model results (default is off to save storage space)
        cfg.boot.write_boot_images = 1; % save out bootstrap thresholded weight maps 
        cfg.boot.write_uncorrected_images = 1; % save out uncorrected statistical maps
        cfg.boot.write_ci_images = 0; % save out statistical maps thresholded by bootstrap CIs (only relevant if bootstrapping CIs)
       
        % Bootstrap P-value options
        cfg.boot.get_pvals = 1; % compute p-values using z-statistics from bootstrap distribution (e.g., Kohoutova et al., 2020 - Nature Protocols)
        cfg.boot.coeff_fwe = 1; % apply family-wise error correction to coefficient p-values (Bonferroni-Holm Method)
        cfg.boot.coeff_fdr = 1; % apply false discovery rate correction to coefficient p-values (Benjamini-Hochberg Method)
        
        % Bootstrap CI options
        cfg.boot.get_cis = 1; % compute confidence intervals on model performance estimates and coefficients using bootstrapped statistics
        cfg.boot.ci_type = 'normal'; % CI estimation method 'normal' for distribution method, 'bcp' for bias-corrected percentile, 'bca' for bias-corrected and accelerated
        cfg.boot.ci_alpha_thresh_model = 0.05; % alpha threshold for CIs on model performance estimates
        cfg.boot.ci_alpha_thresh_coeff = 0.05; % alpha threshold for CIs model coefficients (note: if boot_coeff_fwe == 1, then Bonferroni correction will be applied, but this is only valid for ci_type == 'normal' unless you run enough iterations to allow Bonferroni correction; other CI types based on percentiles are limited by number of resamples)
        cfg.boot.ci_coeff_fwe = 1; % Perform Bonferroni correction on coefficient CIs?
    else
        cfg.boot.write_boot_images = 0;
        cfg.boot.write_uncorrected_images = 0;
    end

    % Permutation options for explanatory model
    if cfg.permutation == 1
        
        % Iterations
        cfg.perm.n_perm = 1000; % number of permutation iterations (1000 is usually good enough unless you want voxel-level p-values; minimum p at 1000 iterations is 0.0009)
        
        if ~contains(cfg.model_spec, {'municorr', 'ttest', 'bmunz', 'munilr', 'muniolsr', 'prop_sub'})
            % Permutation p-value options (cFWE turned on for mass univariate analyses)
            cfg.perm.coeff_cfwe = 0; % Perform continuous FWE correction for voxel permutation p-values, uses value set by fwe_thresh (Mirman et al., 2018 - Neuropsycholgia; note - this is only proven for mass-univariate statistics, but empirically it seems to work well with ridge regression and reasonably well with SVR; for PLS it is observed to be hyper-conservative)
            cfg.perm.coeff_p = 0; % compute voxel-wise permutation p-values (i.e., minimum achievable is 1 / (n_perm + 1))
            cfg.perm.coeff_fdr = 0; % perform FDR correction on voxel-wise permutation p-values (not recommended for voxel-level analyses due to large number of permutations required)
            cfg.perm.coeff_fwe = 0; % perform FWE correction on voxel-wise permutation p-values (not recommended for voxel-level analyses due to large number of permutations required)      
            
            cfg.perm.save_perm_results = 0; % save out permutation model results (default is off to save storage space)
            cfg.perm.write_perm_images = 0; % don't write permutatiom images       
            cfg.perm.write_uncorrected_images = 0; % save out uncorrected statistical maps
        else
            % Permutation p-value options (turned off by defualt for multivariate analyses)
            cfg.perm.coeff_cfwe = 1; % Perform continuous FWE correction for voxel permutation p-values, uses value set by fwe_thresh (Mirman et al., 2018 - Neuropsycholgia; note - this is only proven for mass-univariate statistics, but empirically it seems to work well with ridge regression and reasonably well with SVR; for PLS it is observed to be hyper-conservative)
            cfg.perm.coeff_p = 0; % compute voxel-wise permutation p-values (i.e., minimum achievable is 1 / (n_perm + 1))
            cfg.perm.coeff_fdr = 0; % perform FDR correction on voxel-wise permutation p-values (not recommended for voxel-level analyses due to large number of permutations required)
            cfg.perm.coeff_fwe = 0; % perform FWE correction on voxel-wise permutation p-values (not recommended for voxel-level analyses due to large number of permutations required)      
            
            cfg.perm.save_perm_results = 0; % save out permutation model results (default is off to save storage space)
            cfg.perm.write_perm_images = 1; % don't write permutatiom images       
            cfg.perm.write_uncorrected_images = 0; % save out uncorrected statistical maps
        end
      
    else
        cfg.perm.write_perm_images = 0;
        cfg.perm.write_uncorrected_images = 0;
    end    
   
end

% Within-Dataset Predictive modeling (i.e., nested cross-validation with repeats to assess out-of-sample prediction)
if cfg.cross_validation == 1

    % Cross-validation settings
    cfg.cv.stacked_model = 1; % force inclusion of patients lacking predictor data to enable model stacking
    cfg.cv.cv_type = 'KFold';
    cfg.cv.folds = 5; % number of outer folds for nested cross-validation
    cfg.cv.repeats = 5; % number of repeats for outer loop train/test splits
    cfg.cv.permutation = 0; % permutation testing for nested CV
    cfg.cv.n_perm = 0; % Number of permutation iterations for nested CV permutation test - do minimum necessary, especially if large number of repeats
    cfg.cv.summary_type = 'mean'; % method for summarizing across folds/repeats
    cfg.cv.save_cv_results = 1; % save nested cross-validation results
    cfg.cv.save_perm_cv_results = 0; % save permutation test results (not recommended unless debugging)
    
else
    cfg.cv.stacked_model = 0;
end

% Use cross-validation to optimize hyper-parameters
if cfg.optimize_hyperparams == 1

    % Cross-validation results (for hyper-parameter tuning)
    cfg.hp_opt.cv_type = 'KFold'; % cross-validation, can be KFold or LOOCV
    cfg.hp_opt.folds = 5; % if cv_type is KFold, number of folds 

    switch cfg.model_spec
        case {'linsvr', 'kernsvr', 'linsvc', 'kernsvc'}
            cfg.hp_opt.repeats = 1; % These methods use bayesian optimization with repartitioning at each iteration, repeats are not used
        otherwise
            cfg.hp_opt.repeats = 5; % Other methods use 5 repeats of 5-fold CV by default
    end
end

% P-value thresholds for thresholding and saving out images
cfg.unc_thresh = 0.001; % uncorrected p-value thresholds (only used if write_uncorrected_images == 1)
cfg.fdr_thresh = 0.05; % FDR-corrected p-value thresholds (only used if p-values are bootstrapped and FDR-corrected, and write_boot_images == 1)
cfg.fwe_thresh = 0.05; % FWE-corrected thresholds (used if p-values are bootstrapped and FWE-corrected and write_boot_images == 1)
cfg.direction = 'both'; % Direction for "max" summary
cfg.map_summary = 'max'; % 'max' or 'zscore' - 'max' rescales as a proportion of the maximum (or minimum if cfg.direction == 'neg') value in the map. If cfg.direction == 'both', then the max absolute value is used

%%%% Parallel processing options
% Check for parallel capabilities and turn on if found
my_version = version('-release');
my_version = str2double(my_version(1:end-1));
cfg.parallel = 0;
if my_version >= 2022
    if canUseParallelPool
        cfg.parallel = 1; % use parallel processing
        disp('MATLAB parallel processing functionality detected; parallel processing enabled');
    else
        cfg.parallel = 0; % don't use parallel processing
        disp('MATLAB parallel processing functionality not detected; parallel processing disabled');        
    end
else
    disp('MATLAB version earlier than R2022, checking for Parallel Computing Toolbox');
    available_toolboxes = ver;
    for i = 1:size(available_toolboxes,2)
        if convertCharsToStrings(available_toolboxes(i).Name) == "Parallel Computing Toolbox"
           cfg.parallel = 1;
        end
    end
    if cfg.parallel == 0
     disp('Parallel Computing Toolbox not detected; parallel processing disabled.');
    end
end

% Options to save main / non-CV analysis outputs automatically
cfg.save_model_results = 1; % save out full model results

end
