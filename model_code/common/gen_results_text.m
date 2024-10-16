function result_text = gen_results_text(model_results)

% This function automatically generates a "template" Methods section describing the modeling analysis that was run by the toolkit.
% Joseph Griffis, 2024

% Get cfg as separate structure for convenience
cfg = model_results.cfg;

result_text = [''];

switch cfg.model_spec 
    case {'plsr', 'ridge', 'lasso', 'rlinsvr', 'linsvr', 'kernsvr', 'rensemble'}
        if cfg.fit_explanatory_model == 1
            if cfg.bootstrap == 0
                result_text = ['The inferential model fit to the full dataset explained ' num2str(100*model_results.r2) '% of the variance in the outcome.'];    
            elseif cfg.bootstrap == 1 && cfg.boot.get_cis == 1
                result_text = ['The inferential model fit to the full dataset explained ' num2str(100*model_results.r2) '% of the variance in the outcome ',...
                    '(' num2str(100*(1-cfg.boot.ci_alpha_thresh_model)) '% CI: [' num2str(model_results.boot.r2_ci(1)) ', ' num2str(model_results.boot.r2_ci(2)) ']).'];   
            end
            if cfg.permutation == 1
                if model_results.perm.model_pval <= 0.05
                    result_text = [result_text, ' The permutation test indicated that the model fit better than expected given the permutation null distribution (p=' num2str(model_results.perm.model_pval) ').'];
                else
                    result_text = [result_text, ' The permutation test indicated that the model did not fit better than expected given the permutation null distribution (p=' num2str(model_results.perm.model_pval) ').'];
                end
            end
            if cfg.cross_validation == 1
                load(fullfile(model_results.cfg.out_dir, 'cv_results.mat'));
                if  cv_results.all.corr_pval <= 0.05
                    result_text = [result_text, ' The cross-validation test indicated that predicted out-of-fold outcomes were significantly correlated with the observed outcomes across the full dataset (r=' num2str(cv_results.all.corr) ', p=' num2str(cv_results.all.corr_pval) '), supporting the interpretation of the inferential model.'];
                else
                    result_text = [result_text, ' The cross-validation test indicated that predicted out-of-fold outcomes were not significantly correlated with the observed outcomes across the full dataset (r=' num2str(cv_results.all.corr) ', p=' num2str(cv_results.all.corr_pval) ') and did not provide support for the interpretation of the inferential model.'];
                end
            end
        end
        if cfg.cross_validation == 1
            load(fullfile(model_results.cfg.out_dir, 'cv_results.mat'));            
            result_text = [result_text, ' Across all folds and repeats, the average cross-validation R-squared was ' num2str(mean(cv_results.r2_ss(:))) ' (SD=' num2str(std(cv_results.r2_ss(:))) ').'];
            if cfg.cv.permutation == 1 && cv_results.perm_pval <= 0.05
                result_text = [result_text, ' Permutation testing on the full cross-validation procedure indicated that the average cross-validation loss (MSE) was lower than expected given the permutation null distribution (p=' num2str(cv_results.perm_pval) ').'];
            elseif cfg.cv.permutation == 1 && cv_results.perm_pval > 0.05
                result_text = [result_text, ' Permutation testing on the full cross-validation procedure indicated that the average cross-validation loss (MSE) was not lower than expected given the permutation null distribution (p=' num2str(cv_results.perm_pval) ').'];
            end
        end
    case {'pls_da', 'logistic_ridge', 'logistic_lasso', 'rlinsvc', 'linsvc', 'kernsvc', 'censemble'}
        if cfg.fit_explanatory_model == 1
            if cfg.bootstrap == 0
                result_text = ['The inferential model fit to the full dataset achieved ' num2str(100*model_results.classrate(1)) '% classification accuracy for the group labeled -1, ',...
                    num2str(100*model_results.classrate(2)) '% accuracy for the group labeled 1, and ' num2str(100*model_results.classrate(3)) '% accuracy for the combined groups.',...
                    'The area under the ROC curve for the inferential model fit to the full dataset was ' num2str(model_results.roc_auc) '.'];    
            elseif cfg.bootstrap == 1 && cfg.boot.get_cis == 1
                result_text = ['The inferential model fit to the full dataset achieved ' num2str(100*model_results.classrate(1)) '% classification accuracy for the group labeled -1 ',...
                    '(' num2str(100*(1-cfg.boot.ci_alpha_thresh_model)) '% CI: [' num2str(model_results.boot.classrate_ci(1,1)) ', ' num2str(model_results.boot.classrate_ci(2,1)) ']), ',...
                    num2str(100*model_results.classrate(2)) '% accuracy for the group labeled 1 ',...
                    '(' num2str(100*(1-cfg.boot.ci_alpha_thresh_model)) '% CI: [' num2str(model_results.boot.classrate_ci(1,2)) ', ' num2str(model_results.boot.classrate_ci(2,2)) ']), ',...
                    'and ' num2str(100*model_results.classrate(3)) '% accuracy for the combined groups ',...
                    '(' num2str(100*(1-cfg.boot.ci_alpha_thresh_model)) '% CI: [' num2str(model_results.boot.classrate_ci(1,3)) ', ' num2str(model_results.boot.classrate_ci(2,3)) ']). ',...
                    'The area under the ROC curve for the inferential model fit to the full dataset was ' num2str(model_results.roc_auc) ' (' num2str(100*(1-cfg.boot.ci_alpha_thresh_model)) '% CI: [' num2str(model_results.boot.roc_auc_ci(1)) ', ' num2str(model_results.boot.roc_auc_ci(2)) ']).'];                 
            end                 
            if cfg.permutation == 1
                if model_results.perm.model_pval <= 0.05
                    result_text = [result_text, ' The permutation test indicated that the inferential model fit the full dataset better than expected given the permutation null distribution (p=' num2str(model_results.perm.model_pval) ').'];
                else
                    result_text = [result_text, ' The permutation test indicated that the inferential model did not fit the full dataset better than expected given the permutation null distribution (p=' num2str(model_results.perm.model_pval) ').'];
                end
            end
            if cfg.cross_validation == 1
                load(fullfile(model_results.cfg.out_dir, 'cv_results.mat'));
                if  cv_results.all.fisher_pval <= 0.05
                    result_text = [result_text, ' The cross-validation test indicated that predicted out-of-fold outcomes were significantly associated with the observed outcomes (OR=' num2str(cv_results.all.fisher_stat.OddsRatio) ', p=' num2str(cv_results.all.fisher_pval) '), supporting the interpretation of the inferential model.'];
                else
                    result_text = [result_text, ' The cross-validation test indicated that predicted out-of-fold outcomes were not significantly correlated with the observed outcomes (OR=' num2str(cv_results.all.fisher_stat.OddsRatio) ', p=' num2str(cv_results.all.fisher_pval)  ') and did not support interpretation of the inferential model.'];
                end
            end
        end
        if cfg.cross_validation == 1
            load(fullfile(model_results.cfg.out_dir, 'cv_results.mat'));            
            result_text = [result_text, ' Cross-validation was used to estimate out-of-sample prediction performance. Across all folds and repeats, the average cross-validation area under the ROC curve was ' num2str(mean(cv_results.roc_auc(:))) ' (SD=' num2str(std(cv_results.roc_auc(:))) ').'];
            if cfg.cv.permutation == 1 && cv_results.perm_pval <= 0.05
                result_text = [result_text, ' Permutation testing on the full cross-validation procedure indicated that the average cross-validation area under the ROC curve was greater than expected given the permutation null distribution (p=' num2str(cv_results.perm_pval) ').'];
            elseif cfg.cv.permutation == 1 && cv_results.perm_pval > 0.05
                result_text = [result_text, ' Permutation testing on the full cross-validation procedure indicated that the average cross-validation area under the ROC curve was not greater than expected given the permutation null distribution (p=' num2str(cv_results.perm_pval) ').'];
            end
        end 
    case {'municorr', 'bmunz', 'olsr', 'munilr', 'ttest', 'muniolsr'}
        if isfield(model_results, 'perm')
            if isfield(model_results.perm, 'cfwe')
                crit_val = string(fieldnames(model_results.perm.cfwe));
                if strcmp(model_results.cfg.model_spec, 'bmunz')
                    model_results.tstat = model_results.coeff;
                end
                for i = 1:length(crit_val)
                    n_survive = numel(find(model_results.tstat >= model_results.perm.cfwe.(crit_val(i))));
                    v_thresh = strsplit(crit_val(i), '_');
                    result_text = [result_text, 'For the continuous FWE analyses, the number of surviving voxels at each voxel count threshold (v) was evaluated. At a critical voxel threshold of v=' num2str(v_thresh{end}) ', ' num2str(n_survive) ' voxels survived thresholding. '];
                end
            end
        end                   
end