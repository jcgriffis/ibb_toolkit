function [cv_results, perm_cv_results] = perm_nested_cv(X, Y, cfg_perm)

% Runs permutation tests on nested CV results
% Joseph Griffis 2024

% Define summary function 
if strcmp(cfg_perm.cv.summary_type, 'mean') 
    get_summary = @nanmean;
elseif strcmp(cfg_perm.cv.summary_type, 'median')
    get_summary = @nanmedian;
end

% Preallocate outputs
if cfg_perm.cat_Y == 0    
    mse = zeros(cfg_perm.cv.n_perm,1);

elseif cfg_perm.cat_Y == 1    
    roc_auc = zeros(cfg_perm.cv.n_perm,1);
end

% Run permutation analyses
for i = 1:cfg_perm.cv.n_perm
   
    % Run nested CV using this permutation
    temp_results = run_nested_cv(X, Y, cfg_perm, 1);
    
    % Get results
    if cfg_perm.cat_Y == 0
        mse(i,1) = get_summary(temp_results.avg.mse);
    elseif cfg_perm.cat_Y == 1 
        roc_auc(i) = get_summary(temp_results.avg.roc_auc);
    end
end

% Load original results
cv_results = load(fullfile(cfg_perm.out_dir, 'cv_results.mat'));
cv_results = cv_results.cv_results;

% Compute p-values
if cfg_perm.cat_Y == 0
    mean_mse = get_summary(cv_results.avg.mse(:));
    cv_results.perm_pval = (numel(find(mse < mean_mse))+1) ./ (cfg_perm.cv.n_perm+1);
elseif cfg_perm.cat_Y == 1
    mean_roc_auc = get_summary(cv_results.avg.roc_auc);
    cv_results.perm_pval = (numel(find(roc_auc > mean_roc_auc))+1) ./ (cfg_perm.cv.n_perm+1);
end

% Save p values to CV results
save('cv_results.mat', 'cv_results', '-append');

if isfield(cfg_perm.cv, 'save_perm_cv_results')
    if cfg_perm.cv.save_perm_cv_results == 1
       if cfg_perm.cat_Y == 0
          perm_cv_results.mse = mse;
       else
          perm_cv_results.roc_auc = roc_auc;
       end
       cd(cfg_perm.out_dir);
       if isfile('perm_cv_results.mat')
           save('perm_cv_results.mat',  'perm_cv_results', '-append');
       else
           save('perm_cv_results.mat',  'perm_cv_results', '-v7.3');
       end      
    end
end
end