function output = get_regmdl_hp_opts(cfg, N)

% Hyperparameter Optimization Options
% Joseph Griffis 2024

% Cross-validation Settings
if cfg.hp_opt.bayes_opt == 0
    switch cfg.hp_opt.cv_type
        case 'Holdout'
            if ~isempty(cfg.strat_var) 
                disp('Running repeated Holdout optimization for model hyperparameters');
                disp(['Holdout percentage: ' num2str(cfg.hp_opt.holdout .* 100)]);
                disp(['Number of repeated splits: ' num2str(cfg.hp_opt.repeats)]);
                if isfield(cfg, 'train_set')
                    c = cvpartition(cfg.strat_groups(cfg.train_set), 'Holdout', cfg.hp_opt.holdout);
                else
                    c = cvpartition(cfg.strat_groups, 'Holdout', cfg.hp_opt.holdout);
                end
            else
                disp('Running Holdout optimization for model hyperparameters with repartitioning at each iteration');
                disp(['Holdout percentage: ' num2str(cfg.hp_opt.holdout .* 100)]);
                if isfield(cfg, 'train_set')
                    c = cvpartition(cfg.train_set, 'Holdout', cfg.hp_opt.holdout);
                else
                    c = cvpartition(N, 'Holdout', cfg.hp_opt.holdout);
                end
            end
        case 'KFold' 
            if ~isempty(cfg.strat_var) 
                disp(['Running ' num2str(cfg.hp_opt.repeats) ' repeats of '  num2str(cfg.hp_opt.folds) '-fold optimization for model hyperparameters']);
                if isfield(cfg, 'train_set')
                    c = cvpartition(cfg.strat_groups(cfg.train_set), 'KFold', cfg.hp_opt.folds);
                else
                    c = cvpartition(cfg.strat_groups, 'KFold', cfg.hp_opt.folds);
                end
            else
                disp(['Running ' num2str(cfg.hp_opt.repeats) ' repeats of '  num2str(cfg.hp_opt.folds) '-fold optimization for model hyperparameters']);
                if isfield(cfg, 'train_set')
                    c = cvpartition(cfg.train_set, 'KFold', cfg.hp_opt.folds);
                else
                    c = cvpartition(N, 'KFold', cfg.hp_opt.folds);
                end      
            end
        case 'LOOCV'
            disp('Running leave-one-out optimization for model hyperparameters');            
            c = cvpartition(N, 'LeaveOut');
    end
    output = c;
else
    % Parallel processing
    if cfg.parallel == 1
        hp_opt.UseParallel=true;
    else
        hp_opt.UseParallel=false;
    end
    
    % Defaults for CLI/Plot output
    hp_opt.Verbose=0;
    hp_opt.ShowPlots=false;
    
    % Cross-validation Settings
    switch cfg.hp_opt.cv_type
        case 'Holdout'
            if ~isempty(cfg.strat_var) 
                disp('Running repeated Holdout optimization for model hyperparameters');
                disp(['Holdout percentage: ' num2str(cfg.hp_opt.holdout .* 100)]);
                disp(['Number of repeated splits: ' num2str(cfg.hp_opt.repeats)]);
                if isfield(cfg, 'train_set')
                    c = cvpartition(cfg.strat_groups(cfg.train_set), 'Holdout', cfg.hp_opt.holdout);
                else
                    c = cvpartition(cfg.strat_groups, 'Holdout', cfg.hp_opt.holdout);
                end
                hp_opt.CVPartition=c;
            else
                disp('Running Holdout optimization for model hyperparameters with repartitioning at each iteration');
                disp(['Holdout percentage: ' num2str(cfg.hp_opt.holdout .* 100)]);
                hp_opt.Holdout=cfg.hp_opt.holdout;
                hp_opt.NumGridDivisions=cfg.hp_opt.grid_div;
                hp_opt.Repartition=true;
                hp_opt.Optimizer='gridsearch';
            end
        case 'KFold' 
            if ~isempty(cfg.strat_var) 
                disp('Running repeated KFold optimization for model hyperparameters');
                disp(['Number of folds: ' num2str(cfg.hp_opt.folds)]);
                disp(['Number of repeated ' num2str(cfg.hp_opt.folds) '-fold splits: ' num2str(cfg.hp_opt.repeats)]);
                if isfield(cfg, 'train_set')
                    c = cvpartition(cfg.strat_groups(cfg.train_set), 'KFold', cfg.hp_opt.folds);
                else
                    c = cvpartition(cfg.strat_groups, 'KFold', cfg.hp_opt.folds);
                end
                hp_opt.CVPartition=c;
            else
                disp('Running KFold optimization for model hyperparameters');
                disp(['Number of folds: ' num2str(cfg.hp_opt.folds)]);
                disp(['Number of grid divisions: ' num2str(cfg.hp_opt.grid_div)]);   
                if cfg.hp_opt.repartition == true
                    disp('Repartitioning at each iteration');
                end
                hp_opt.Kfold=cfg.hp_opt.folds;
                hp_opt.Repartition=cfg.hp_opt.repartition;
                hp_opt.NumGridDivisions = cfg.hp_opt.grid_div;
                hp_opt.Optimizer = 'gridsearch';
            end
        case 'LOOCV'
            c = cvpartition(N, 'LeaveOut');
            hp_opt.CVPartition=c;
    end
    
    % Acquisition function name
    if ~isfield(cfg, 'AcqFuncName')
        hp_opt.AcquisitionFunctionName='expected-improvement-plus';
    else
        hp_opt.AcquisitionFunctionName=cfg.hp_opt.acq_func_name;
    end
    output = hp_opt;
end
end