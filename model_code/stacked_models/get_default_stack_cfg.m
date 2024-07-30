function cfg = get_default_stack_cfg(cfg)

% Generates a pre-filled cfg structure containing modeling parameters and analysis options
% Joseph Griffis 2024

% Add common directory to path
restore_model_path(cfg.model_dir);
addpath(fullfile(cfg.model_dir, 'common'));

% Set defaults if not already set
if ~isfield(cfg, 'model_spec')
   error('cfg.model_spec is not set. Please specify the model and try again.')
end
if ~isfield(cfg, 'optimize_hyperparams')
    disp('cfg.optimize_hyperparams field is empty; setting to 1 with default settings')
    cfg.optimize_hyperparams = 1;
end
if ~isfield(cfg, 'cross_validation')
    disp('cfg.cross_validation field is empty; setting to 1 with default settings')
    cfg.cross_validation = 1;
end

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
        cfg.use_scores = 1;        
    case {'lasso', 'ridge', 'rlinsvr', 'logistic_lasso', 'logistic_ridge', 'rlinsvc'} % Regularized linear models
        if strcmp(cfg.model_spec, {'lasso', 'ridge'})
            cfg.learner = 'leastsquares'; % OLS regression
        elseif contains(cfg.model_spec, 'logistic')
            cfg.learner = 'logistic'; % logistic regression
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
        cfg.use_scores = 1;        
        addpath(fullfile(cfg.model_dir, 'reglm'));              
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
        cfg.use_scores = 1;        
        cfg.hp_opt.bayes_opt = 1; % Do Bayesian hyper-parameter optimization
        cfg.hp_opt.opt_iter = 100; % Number of objective function evaluations
        cfg.hp_opt.repartition = true; % Repartition into train/test at each iteration 
        addpath(fullfile(cfg.model_dir, 'nlinsvr'));
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
        addpath(fullfile(cfg.model_dir, 'ensemble'));         
    case 'mean'
        cfg.cat_Y = 0;
        cfg.use_scores = 1;
        addpath(fullfile(cfg.model_dir, 'stacked_models'))
    case 'mean_score' 
        cfg.cat_Y = 1;
        cfg.use_scores = 1;
        addpath(fullfile(cfg.model_dir, 'stacked_models'))
    case 'olsr'
        cfg.cat_Y = 0;
        addpath(fullfile(cfg.model_dir, 'mass_univariate'));
    case 'munilr' % Mass univariate logistic regression - not currently functional 
        cfg.cat_Y = 1;
        cfg.use_scores = 1;
        cfg.hp_opt = [];        
        addpath(fullfile(cfg.model_dir, 'mass_univariate'));                     
    otherwise
        warning('Model type not provided as input - make sure to set model-specific CFG fields since they will not be set by default!' );
end
if cfg.cat_Y == 1
    cfg.cost = [0,1;1,0];
end

% Confounds
cfg.confounds = [];

% Standardization
cfg.standardize = 0;

% Stratification
cfg.strat_var = [];

% Use cross-validation to optimize hyper-parameters
if cfg.optimize_hyperparams == 1

    % Cross-validation results (for hyper-parameter tuning)
    cfg.hp_opt.cv_type = 'KFold'; % cross-validation, can be KFold or LOOCV
    cfg.hp_opt.folds = 5; % if cv_type is KFold, number of folds 

    if strcmp(cfg.model_spec, 'plsr')
        cfg.hp_opt.repeats = 5; % repeats for CV
    else
        cfg.hp_opt.repeats = 1;
    end
end

% Save results
cfg.cv.save_cv_results = 1;

% Permutation testing
cfg.cv.permutation = 0;
cfg.cv.n_perm = 0;

%%%% Parallel processing options
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

end
