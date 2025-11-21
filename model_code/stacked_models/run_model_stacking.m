function [cv_results, base_results] = run_model_stacking(cfg)

% Make output directory if it doesn't exist
if ~isfolder(cfg.out_dir)
    mkdir(cfg.out_dir);
end

% Check whether to run regression or classification models
if cfg.cat_Y == 1
    [cv_results, base_results] = stack_class_models_and_run_cv(cfg);
elseif cfg.cat_Y == 0
    [cv_results, base_results] = stack_reg_models_and_run_cv(cfg);
end

% Do permutation testing if indicated
if cfg.cv.permutation == 1
    cfg.perm_flag = 1; 
    cv_results = perm_stacked_cv(base_results, cv_results);
end

% Save CV results
if cfg.cv.save_cv_results == 1
    cd(cfg.out_dir);
    if isfile('cv_results.mat')
        save('cv_results.mat', 'cv_results', 'base_results', '-append');
    else
        save('cv_results.mat', 'cv_results', 'base_results', '-v7.3');
    end
end

if cfg.fit_explanatory_model == 1

    for i = 1:length(cfg.model_paths)

        % Load results and create full-sample X
        load(fullfile(cfg.model_paths(i), 'cv_results.mat'));
        load(fullfile(cfg.model_paths(i), 'model_results.mat'));

        if cfg.cat_Y == 0
            cfg.X(:,i) = nanmean(cv_results.all.pred_y,2);
        elseif cfg.cat_Y == 1
            if isfield(cfg, 'use_scores')
                cfg.X(:,i) = nanmean(cv_results.all.pred_score,2);
            elseif isfield(cfg, 'use_labels')
                cfg.X(:,i) = mode(cv_results.all.pred_y,2);
            else 
                cfg.X(:,i) = mode(cv_results.all.pred_y,2);
            end
        end
    end

    % Add other predictors as needed
    if isfield(cfg, 'other_predictiors')
        cfg.X = [cfg.X, cfg.other_predictors];
    end

    % Get outcome 
    cfg.Y = model_results.obs_y;

    % Set params to pass trimming 
    cfg.freq_thresh = 0;
    cfg.min_obs = 0;
    cfg.trim_X = 0;
    cfg.gen_freq_map = 0;
    cfg.gen_lvol_map = 0;
    cfg.dtlvc = 0;
    cfg.strat_var = cfg.Y;
    cfg.cross_validation = 0;
    if isfield(model_results.cfg, 'strat_var')
        cfg.strat_var = model_results.cfg.strat_var;
        if isfield(model_results.cfg, 'strat_groups')
            cfg.strat_groups = model_results.cfg.strat_groups;
        end
    end
    clear model_results cv_results

    % Run full-sample model
    cfg.save_model_results = 1;
    fit_and_evaluate_model(cfg);
    
end 

end
