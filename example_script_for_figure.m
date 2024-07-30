% Add paths to relevant functions
cfg.model_dir = 'path/to/repo/mlsm_code_repo/model_code';
addpath(genpath(fullfile(cfg.model_dir, 'common'))); 

% Specify data modality and supply brain mask
cfg.modality = 'lesion'; 
cfg.mask_path = fullfile('path/to/brain_mask.nii'); % path to brain mask

% Specify general modeling options
cfg.model_spec = 'plsr'; 
cfg.fit_explanatory_model = 1; 
cfg.optimize_hyperparams = 1; 
cfg.cross_validation = 1; 
cfg.bootstrap = 1; 
cfg.permutation = 1; 

% Initialize cfg file with default settings given general modeling options
cfg = get_default_model_cfg(cfg); 

% Specify path to CSV containing study IDs and behavioral scores
beh_csv = readtable('path/to/behavior/test_data.csv');
id_col = 'study_id';
beh_col = 'mae_token';

% Specify lesion directory and indicate whether patients are registry patients
lesion_dir = '/path/to/lesion/Lesion_3mm'; % directory containing lesion images
registry_flag = 1; 

% Get formatted lesion and behavioral daa
cfg = get_and_format_lesion_data(beh_csv, id_col, beh_col, registry_flag, lesion_dir, cfg);

% Variable transformation
cfg.dtlvc = 1; 
cfg.standardize = 0;

% Statistical significance testing 
cfg.boot.n_boot = 1000; 
cfg.perm.n_perm = 1000;

% Hyper-parameter optimization
cfg.hp_opt.cv_type = 'KFold'; 
cfg.hp_opt.folds = 5; 
cfg.hp_opt.repeats = 5;

% Cross-validation
cfg.cv.cv_type = 'KFold'; 
cfg.cv.folds = 5; 
cfg.cv.repeats = 5; 
cfg.cv.permutation = 1; 
cfg.cv.n_perm = 50; 
cfg.cv.stacked_model = 1; 

% Stratification for cross-validation 
cfg.strat_var = cfg.Y; % stratification variable

% Specify output directory
cfg.out_dir = fullfile('/path/to/output', [cfg.model_spec '_' beh_col]);

% fit model with CV hyperparameter optimization and selected significance testing
[model_results] = fit_and_evaluate_model(cfg);