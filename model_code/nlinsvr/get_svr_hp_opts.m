function hp_opt = get_svr_hp_opts(cfg, N)

% Hyperparameter Optimization Options
% Joseph Griffis 2024

% Parallel processing
if cfg.parallel == 1
    hp_opt.UseParallel=true;
else
    hp_opt.UseParallel=false;
end

% Defaults for CLI/Plot output
hp_opt.Verbose=1;
hp_opt.ShowPlots=false;

% Cross-validation Settings
switch cfg.hp_opt.cv_type
    case 'KFold' 
        if ~isempty(cfg.strat_var) && cfg.hp_opt.repartition == 0
            if isfield(cfg, 'train_set')
                c = cvpartition(cfg.strat_groups(cfg.train_set), 'KFold', cfg.hp_opt.folds);
            else
                c = cvpartition(cfg.strat_groups, 'KFold', cfg.hp_opt.folds);
            end
            hp_opt.CVPartition=c;
            hp_opt.MaxObjectiveEvaluations=cfg.hp_opt.opt_iter;            
        else
            if cfg.hp_opt.repartition == true
                disp('Repartitioning at each iteration');
            end
            hp_opt.Kfold=cfg.hp_opt.folds;
            hp_opt.Repartition=cfg.hp_opt.repartition;
            hp_opt.MaxObjectiveEvaluations=cfg.hp_opt.opt_iter;
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

end