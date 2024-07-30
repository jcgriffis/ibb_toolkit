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

end