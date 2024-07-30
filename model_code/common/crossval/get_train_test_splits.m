function [train_ind, test_ind] = get_train_test_splits(Y,cfg)

% Generates train and test splits based on options specified in cfg
% Y: response vector (Nx1)
% cfg: structure containing parameters/options
% Joseph Griffis 2023

% Stratify if indicated and partition data
if isfield(cfg, 'strat_var') && ~isempty(cfg.strat_var)
    groups = get_quartile_groups(cfg.strat_var);
    train_ind = zeros(length(Y), cfg.n_repeats, cfg.folds);
    test_ind = zeros(length(Y), cfg.n_repeats, cfg.folds);
    for i = 1:length(cfg.n_repeats)
        c = cvpartition(groups, 'KFold', cfg.folds);
        for j = 1:length(cfg.folds)
            train_ind(:,i,j) = training(c,j);
            test_ind(:,i,j) = test(c,j);
        end
    end

else
    train_ind = zeros(length(Y), cfg.n_repeats, cfg.folds);
    test_ind = zeros(length(Y), cfg.n_repeats, cfg.folds);
    for i = 1:length(cfg.n_repeats)
        c = cvpartition(length(Y), 'KFold', cfg.folds);
        for j = 1:length(cfg.folds)
            train_ind(:,i,j) = training(c,j);
            test_ind(:,i,j) = test(c,j);
        end
    end
end
end
