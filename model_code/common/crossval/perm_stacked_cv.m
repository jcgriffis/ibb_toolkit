function [cv_results, perm_cv_results] = perm_stacked_cv(base_results, cv_results)

% Runs permutation tests on nested CV results
% Joseph Griffis 2024

% Define summary function 
if strcmp(cv_results.cfg.cv.summary_type, 'mean') 
    get_summary = @nanmean;
elseif strcmp(cv_results.cfg.cv.summary_type, 'median')
    get_summary = @nanmedian;
end

% Preallocate outputs
if cv_results.cfg.cat_Y == 0    
    mse = zeros(cv_results.cfg.cv.n_perm,1);
elseif cv_results.cfg.cat_Y == 1    
    roc_auc = zeros(cv_results.cfg.cv.n_perm,1);
end

% Run permutation analyses
for i = 1:cv_results.cfg.cv.n_perm
   
    % Preallocate outputs
    temp_results = preallocate_cv_results_perm(length(base_results.all.pred_y), cv_results.cfg);
    
    % Get permuted Y 
    perm_Y = cv_results.all.obs_y(randperm(length(cv_results.all.obs_y))');

    % Run nested CV using this permutation
    temp_results = run_stacked_cv(perm_Y, cv_results.cfg, temp_results, base_results, 1);
    
    % Get results
    if cv_results.cfg.cat_Y == 0
        mse(i,1) = get_summary(temp_results.avg.mse);
    elseif cv_results.cfg.cat_Y == 1 
        roc_auc(i) = get_summary(temp_results.avg.roc_auc);
    end
end

% Compute p-values
if cv_results.cfg.cat_Y == 0
    mean_mse = get_summary(cv_results.avg.mse(:));
    cv_results.perm_pval = (numel(find(mse < mean_mse))+1) ./ (cv_results.cfg.cv.n_perm+1);
elseif cv_results.cfg.cat_Y == 1
    mean_roc_auc = get_summary(cv_results.avg.roc_auc);
    cv_results.perm_pval = (numel(find(roc_auc > mean_roc_auc))+1) ./ (cv_results.cfg.cv.n_perm+1);
end

% Save p values to CV results
save('cv_results.mat', 'cv_results', '-append');

if isfield(cv_results.cfg.cv, 'save_perm_cv_results')
    if cv_results.cfg.cv.save_perm_cv_results == 1
       if cfg.cat_Y == 0
          perm_cv_results.mse = mse;
       else
          perm_cv_results.classrate1 = classrate1;
          perm_cv_results.classrate2 = classrate2;
          perm_cv_results.roc_auc = roc_auc;
       end
       cd(cv_results.cfg.out_dir);
       if isfile('perm_cv_results.mat')
           save('perm_cv_results.mat',  'perm_cv_results', '-append');
       else
           save('perm_cv_results.mat',  'perm_cv_results', '-v7.3');
       end      
    end
end
end