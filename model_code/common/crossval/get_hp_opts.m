function output = get_hp_opts(cfg, N)

% Hyperparameter Optimization Options
% Joseph Griffis 2024

% Cross-validation Settings
if cfg.hp_opt.bayes_opt == 0
    switch cfg.hp_opt.cv_type
        case 'Holdout'
            if ~isempty(cfg.strat_var)
                disp(['Running ' num2str(cfg.hp_opt.repeats) ' repeats of ' num2str(cfg.hp_opt.holdout .* 100) '% hold-out optimization for model hyperparameters']);                
                if isfield(cfg, 'train_set')
                    c = cvpartition(cfg.strat_groups(cfg.train_set==1), 'Holdout', cfg.hp_opt.holdout);
                else
                    c = cvpartition(cfg.strat_groups, 'Holdout', cfg.hp_opt.holdout);
                end
            else
                disp(['Running ' num2str(cfg.hp_opt.repeats) ' repeats of ' num2str(cfg.hp_opt.holdout .* 100) '% hold-out optimization for model hyperparameters']);                
                c = cvpartition(N, 'KFold', cfg.hp_opt.folds);
            end
        case 'KFold' 
            if ~isempty(cfg.strat_var) 
                disp(['Running ' num2str(cfg.hp_opt.repeats) ' repeats of '  num2str(cfg.hp_opt.folds) '-fold optimization for model hyperparameters']);
                if isfield(cfg, 'train_set')
                    c = cvpartition(cfg.strat_groups(cfg.train_set==1), 'KFold', cfg.hp_opt.folds);
                else
                    c = cvpartition(cfg.strat_groups, 'KFold', cfg.hp_opt.folds);
                end
            else
                disp(['Running ' num2str(cfg.hp_opt.repeats) ' repeats of '  num2str(cfg.hp_opt.folds) '-fold optimization for model hyperparameters']);
                c = cvpartition(N, 'KFold', cfg.hp_opt.folds);
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
            disp('Running Holdout optimization for model hyperparameters with repartitioning at each iteration');
            disp(['Holdout percentage: ' num2str(cfg.hp_opt.holdout .* 100)]);
            disp(['Number of optimization iterations: ' num2str(cfg.hp_opt.opt_iter)]);   
            hp_opt.Holdout=cfg.hp_opt.holdout;
            hp_opt.MaxObjectiveEvaluations=cfg.hp_opt.opt_iter;
            if cfg.hp_opt.repartition == true
                disp('Repartitioning at each iteration');
            end
            hp_opt.Repartition=cfg.hp_opt.repartition;            
        case 'KFold' 
            disp('Running KFold optimization for model hyperparameters');
            disp(['Number of folds: ' num2str(cfg.hp_opt.folds)]);
            disp(['Number of optimization iterations: ' num2str(cfg.hp_opt.opt_iter)]);   
            if cfg.hp_opt.repartition == true
                disp('Repartitioning at each iteration');
            end
            hp_opt.Kfold=cfg.hp_opt.folds;
            hp_opt.Repartition=cfg.hp_opt.repartition;
            hp_opt.MaxObjectiveEvaluations=cfg.hp_opt.opt_iter;
        case 'LOOCV'
            disp('Running leave-one-out optimization for model hyperparameters');
            disp(['Number of optimization iterations: ' num2str(cfg.hp_opt.opt_iter)]);               
            c = cvpartition(N, 'LeaveOut');
            hp_opt.CVPartition=c;
            hp_opt.MaxObjectiveEvaluations=cfg.hp_opt.opt_iter;            
    end
    
    % Acquisition function name
    if isfield(cfg, 'AcqFuncName')
       hp_opt.AcquisitionFunctionName=cfg.hp_opt.acq_func_name;
    end
    output = hp_opt;
end
end
