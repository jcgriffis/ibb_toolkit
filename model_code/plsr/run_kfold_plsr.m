function [opt_k, all_mse] = run_kfold_plsr(X,Y,cfg)
% This function identifies the optimal number of PLSR components
% Inputs:
% X: predictor matrix (Nxk)
% Y: response vector (Nx1)
% cfg: structure containing parameters/options
% Joseph Griffis 2023

% Max number of components
if isfield(cfg.hp_opt, 'n_comp')
    disp(['Using pre-specified upper bound on components: ' num2str(cfg.hp_opt.n_comp) ' components']);
    if cfg.hp_opt.n_comp <= (min(size(X))-1)
        n_comp = cfg.hp_opt.n_comp;
    else
        n_comp = min(size(X))-1;
    end
else
    disp(['Using maximum possible number of: ' num2str(cfg.hp_opt.n_comp) ' components']);
    n_comp = min(size(X))-1;
end    

% Parallel flag
if cfg.parallel == 1
    options = statset('UseParallel',true);
else
    options = statset('UseParallel',false);
end

% Check to make sure repeats is non-zero
if isfield(cfg.hp_opt, 'repeats') && cfg.hp_opt.repeats == 0
    cfg.hp_opt.repeats = 1;
elseif ~isfield(cfg.hp_opt, 'repeats')
    cfg.hp_opt.repeats = 1;
end

%%%% Repeat CV and estimate optimal component number for each CV

% For stratified sampling, can't use PLSREGRESS internal CV
if ~isempty(cfg.strat_var)

    % Stratify if indicated and partition data
    if isfield(cfg, 'train_set') && ~isempty(cfg.train_set) % if this is nested CV
        c = cvpartition(cfg.strat_groups(cfg.train_set), 'KFold', cfg.hp_opt.folds, 'Stratify', true);
    else
        c = cvpartition(cfg.strat_groups, 'KFold', cfg.hp_opt.folds, 'Stratify', true);
    end

    % pre-allocate outputs
    all_mse = zeros(length(cfg.hp_opt.repeats), n_comp);

    for i = 1:cfg.hp_opt.repeats
        disp(['Running cross-validation iteration: ' num2str(i) '/' num2str(cfg.hp_opt.repeats)]);
        % repartition data 
        if i > 1 
            c = c.repartition;
        end
        
        % Make sure n_comp is valid given size of training set
        for j = 1:length(c.NumTestSets)
            n_train = numel(find(training(c,j)));
            if n_comp > (n_train - 1)
                n_comp = n_train - 1;
            else
                n_comp = cfg.hp_opt.n_comp;
            end
        end

        % Fit PLS model with Kfold cross-validation
        [~,~,~,~,~,~, MSE] = plsregress(X, Y, n_comp,'cv', c, 'Options', options);
        
        % Get MSE for outcome
        y_mse = MSE(2,2:end);

        if length(y_mse) < n_comp
            % Buffer with NaNs if sample is too small to get n_comp components
            len_dif = n_comp - length(y_mse);
            y_mse((end-len_dif):end) = NaN;
        end

        all_mse(i,:) = y_mse;

    end

    % Get average MSE across repeats
    if cfg.hp_opt.repeats > 1
        mean_mse = mean(all_mse);
    else
        mean_mse = all_mse;
    end
    
    switch cfg.hp_opt.comp_method
       % Find index of last component that decreased cv-MSE
       case 'dif_mse'
           diff_mse = diff(mean_mse);
           mse_reversal = find(diff_mse > 0, 1);        
           if ~isempty(mse_reversal)
              opt_k = mse_reversal;
           else
              opt_k = n_comp;
           end
       % Find index of component with minimum CV-MSE out of all n_comp
       case 'min_mse'
           [~, min_mse] = min(mean_mse);
           opt_k = min_mse(1);
    end
    if cfg.hp_opt.repeats > 1
        all_mse = mean_mse;
    end

% If no stratification, use PLSREGRESS internal CV 
else

    disp(['Using PLSRegress internal CV with:' num2str(cfg.hp_opt.repeats) ' repeats']);
    % Fit PLS model with specified folds and repeats
    [~,~,~,~,~,~, MSE] = plsregress(X, Y, n_comp, 'cv', cfg.hp_opt.folds, 'mcreps', cfg.hp_opt.repeats, 'Options',options);
    
   switch cfg.hp_opt.comp_method
       % Find index of last component that decreased cv-MSE
       case 'dif_mse'
           all_mse = MSE(2,2:end);
           diff_mse = diff(all_mse);
           mse_reversal = find(diff_mse > 0, 1);        
           if ~isempty(mse_reversal)
              opt_k = mse_reversal;
           else
              opt_k = n_comp;
           end
       % Find index of component with minimum CV-MSE out of all n_comp
       case 'min_mse'
           all_mse = MSE(2,2:end);
           [~, min_mse] = min(all_mse);
           opt_k = min_mse(1);
    end
end

end

