function [stacked_results, base_results] = stack_class_models_and_run_cv(cfg)

% This function takes existing cross-validation predictions from different models, 
% and re-runs the cross-validation to evaluate a stacked model combining the predictions 
% across models into a single model

% Joseph Griffis 2024

% Go to model paths and load base model results
base_results.all.pred_y = cell(size(cfg.model_paths));
base_results.all.pred_score = cell(size(cfg.model_paths));
base_results.avg.classrate = cell(size(cfg.model_paths));
base_results.avg.roc_auc = cell(size(cfg.model_paths));
base_results.pred_y = cell(size(cfg.model_paths));
base_results.classrate = cell(size(cfg.model_paths));
base_results.roc_auc = cell(size(cfg.model_paths));

for i = 1:length(cfg.model_paths)

    % Load results
    load(fullfile(cfg.model_paths(i), 'cv_results.mat'));
    if i == 1
        load(fullfile(cfg.model_paths(i), 'model_results.mat'));
    end

    % Get performance measures collapsed across test folds
    base_results.avg.classrate{i} = cv_results.avg.classrate;
    base_results.avg.roc_auc{i} = cv_results.avg.roc_auc;
    base_results.all.pred_y{i} = cv_results.all.pred_y;    
    base_results.all.pred_score{i} = cv_results.all.pred_score;

    % Get fold-level performance measures 
    base_results.classrate{i} = cv_results.classrate;
    base_results.roc_auc{i} = cv_results.roc_auc;
    base_results.pred_y{i} = cv_results.pred_y;

    % Get relevant data    
    if i == 1
        cfg.cross_validation = 1;
        cfg.cv.cv_type = model_results.cfg.cv.cv_type;
        cfg.cv.repeats = size(cv_results.train_set,2);
        cfg.cv.folds = size(cv_results.train_set, 3);        
        cfg.cv.partitions.train_set = cv_results.train_set;
        cfg.cv.partitions.test_set = cv_results.test_set;    
        if isfield(model_results.cfg.cv, 'summary_type')
           cfg.cv.summary_type = model_results.cfg.cv.summary_type;
        else
           cfg.cv.summary_type = 'mean';
        end        
        Y = cv_results.all.obs_y(:,1);
        clear model_results
    else
        assert(isequal(cfg.cv.partitions.train_set, cv_results.train_set), 'Error: Training datasets must be the same to allow for model stacking.')
    end

end 

% Preallocate stacked results output
stacked_results = preallocate_cv_results(length(Y), length(base_results.pred_y), cfg);

% Run cross-validation to fit stacked models
perm_flag = 0;
stacked_results = run_stacked_cv(Y, cfg, stacked_results, base_results, perm_flag);
stacked_results.cfg = cfg;    

end