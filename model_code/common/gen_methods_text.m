function method_text = gen_methods_text(model_results)

% This function automatically generates a "template" Methods section describing the modeling analysis that was run by the toolkit.
% Joseph Griffis, 2024
% (modified by Stein Acker)

% Get cfg as separate structure for convenience
cfg = model_results.cfg;

if ~isfield(cfg, 'strat_name') && cfg.strat_var ~= model_results.obs_y
    cfg.strat_name = '[INSERT VARIABLE NAME]';
elseif cfg.strat_var == model_results.obs_y
    cfg.strat_name = ' the outcome';
end

if ~isfield(cfg.hp_opt, 'bayes_opt')
    cfg.hp_opt.bayes_opt = 0;
end

% Modeling method
switch cfg.model_spec
    case 'plsr'
        mdl_text = 'Data were modeled using partial least squares regression (PLSR).';
    case 'pls_da'
        mdl_text = 'Data were modeled using partial least squares discriminant analysis (PLS-DA).';
    case 'ridge'
        mdl_text = 'Data were modeled using L2-regularized (ridge) regression.';
    case 'lasso'
        mdl_text = 'Data were modeled using L1-regularized (LASSO) regression.';
    case 'rlinsvr'
        mdl_text = 'Data were modeled using L2-regularized (ridge) linear support vector regression (SVR).';
    case 'linsvr'
        mdl_text = 'Data were modeled using linear support vector regression (SVR).'; 
    case 'linsvc'
        mdl_text = 'Data were modeled using L2-regularized (ridge) linear support vector classifier (SVC).';        
    case 'logistic_ridge'
        mdl_text = 'Data were modeled using L2-regularized (ridge) logistic regression.';
    case 'logistic_lasso'
        mdl_text = 'Data were modeled using L1-regularized (LASSO) logistic regression.';
    case 'rlinsvc'
        mdl_text = 'Data were modeled using an L2-regularized (ridge) linear support vector classifier (SVC).';
    case 'kernsvr'
        mdl_text = 'Data were modeled using a support vector regression (SVR) model with a radial basis function (RBF) kernel.';
    case 'kernsvc'
        mdl_text = 'Data were modeled using a support vector classification (SVC) model with a radial basis function (RBF) kernel.';
    case 'municorr'
        mdl_text = 'Data were analyzed using mass-univariate Pearson correlation analyses.';
    case 'ttest'
        mdl_text = 'Data were analyzed using mass-univariate independent samples t-tests.';
    case 'munilr'
        mdl_text = 'Data were analyzed using mass-univariate logistic regression models.';
    case 'bmunz'
        mdl_text = 'Data were analyzed using mass-univariate Brunner-Munzel tests.';
    case 'censemble'
        mdl_text = 'Data were modeled using a classification ensemble.';
    case 'rensemble'
        mdl_text = 'Data were modeled using a regression ensemble';

end

% Analysis configuation
if cfg.cv.stacked_model == 0
    mdl_text = [mdl_text ' Only predictors with at least ' num2str(cfg.freq_thresh) ' non-zero observations with absolute values greater than ' num2str(cfg.min_obs) ' were included in the analysis.'...
        ' After excluding observations with no data for any included predictors, N=' num2str(length(model_results.obs_y)) ' observations were included in the analysis.'];
else
    mdl_text = [mdl_text ' Only predictors with at least ' num2str(cfg.freq_thresh) ' non-zero observations with absolute values greater than ' num2str(cfg.min_obs) ' were included in the analysis.'...
        ' N=' num2str(length(model_results.obs_y)) ' observations were included in the analysis.'];
end
if cfg.dtlvc == 1
    mdl_text = [mdl_text, ' The direct total lesion volume control method (Zhang et al., 2014) was used to mitigate the effects of lesion volume.'];
end
if ~isempty(cfg.confounds)
    if ~isfield(cfg, 'confound_names')
        mdl_text = [mdl_text, ' Confound regression was performed to remove variance associated with user-defined nuisance regressors from the outcome variable in the training set(s) prior to training the model, and the resulting model(s) were applied to the outcome data from the test set(s) before obtaining predicted outcomes.'];
    else
        for i = 1:length(cfg.confound_names)
            if i == 1
                conf_text = cfg.confound_names{i};
            else
                conf_text = [conf_text ', ' cfg.confound_names{i}];
            end
        end
        mdl_text = [mdl_text, ' Confound regression was performed to remove variance associated with user-defined nuisance regresors (' conf_text ') from the outcome variable in the training set(s) prior to training the model, and the resulting model(s) were applied to the outcome data from the test set(s) before obtaining predicted outcomes.'];
    end
end
if cfg.standardize == 1
    mdl_text = [mdl_text, ' Predictor variables were standardized to have means equal to 0 and standard deviations equal to 1.'];
elseif cfg.standardize == 2
    mdl_text = [mdl_text, ' Outcome variables were standardized to have means equal to 0 and standard deviations equal to 1.'];
elseif cfg.standardize == 3
    mdl_text = [mdl_text, ' Predictor and outcome variables were standardized to have means equal to 0 and standard deviations equal to 1.'];    
end


% Cross-validation methods
if cfg.cross_validation == 1
    if cfg.optimize_hyperparams == 1
        cv_text = ' Nested cross-validation was used to optimize model hyper-parameters and obtain estimates of out-of-sample prediction performance.';
    else
        cv_text = ' Cross-validation was used to obtain estimates of out-of-sample prediction performance.';
    end
    if strcmp(cfg.cv.cv_type, 'KFold')
        if cfg.cv.repeats > 1 
            cv_text = [cv_text ' Cross-validation analyses were performed using ' num2str(cfg.cv.repeats) ' repetitions of '...
                num2str(cfg.cv.folds) '-fold cross-validation with repartitioning at each repetition.'];
            if ~isempty(cfg.strat_var)
                cv_text = [cv_text, ' Train and test partitions were stratified according to ' cfg.strat_name '.'];
            end
        else
            cv_text = [cv_text ' Cross-validation analyses were performed using '...
                num2str(cfg.cv.folds) '-fold cross-validation with repartitioning at each repetition.'];
            if ~isempty(cfg.strat_var)
                cv_text = [cv_text, ' Train and test partitions were stratified according to ' cfg.strat_name '.'];
            end
        end
    elseif strcmp(cfg.cv.cv_type, 'Holdout')
        if cfg.cv.repeats > 1 
            cv_text = [cv_text ' Cross-validation analyses were performed using ' num2str(cfg.cv.repeats) ' repetitions of '...
                num2str(cfg.cv.holdout) '% holdout cross-validation with repartitioning at each repetition.'];
            if ~isempty(cfg.strat_var)
                cv_text = [cv_text, ' Train and test partitions were stratified according to ' cfg.strat_name '.'];
            end
        else
            cv_text = [cv_text ' Cross-validation analyses were performed using '...
                num2str(cfg.hp_opt.holdout .* 100) 'holdout cross-validation with repartitioning at each repetition.'];
            if ~isempty(cfg.strat_var)
                cv_text = [cv_text, ' Train and test partitions were stratified according to ' cfg.strat_name '.'];
            end
        end        
    elseif strcmp(cfg.cv.cv_type, 'LOOCV')
        cv_text = [cv_text ' Cross-validation analyses were performed using leave-one-out cross-validation.'];
    end
    if cfg.fit_explanatory_model == 1
        if cfg.cat_Y == 1
            if cfg.cv.repeats > 1
                cv_text = [cv_text, ' To test the significance of the inferential model in terms of the cross-validation performance, Fishers exact test was performed on the confusion matrix computed using the '...
                    'mode (across cross-validation repetitions) out-of-sample predicted classes and true outcomes for all observations (analogous to Pustina et al., 2017).'];
            else
                cv_text = [cv_text, ' To test the significance of the inferential model in terms of cross-validation performance, Fishers exact test was performed on the confusion matrix computed using the '...
                    'out-of-fold predicted classes and true outcomes for all observations (analogous to Pustina et al., 2017).'];            
            end
        else
            if cfg.cv.repeats > 1
                cv_text = [cv_text, ' To test the significance of the inferential model in terms of the cross-validation performance, the average (across cross-validation repetitions) Pearson correlation was computed between the out-of-fold predictions and observed outcomes for all observations (Pustina et al., 2017).'];
            else
                cv_text = [cv_text, ' To test the significance of the inferential model in terms of the cross-validation performance, the Pearson correlation was computed between the out-of-fold predictions and the observed outcomes for all observations (Pustina et al., 2017).'];
            end
        end
    end       
    if cfg.cat_Y == 0
        cv_text = [cv_text ' Out-of-sample prediction performance was estimated using the cross-validation explained variance score, prediction R-squared, Pearson correlation of predicted and observed outcomes, and mean squared error. '];
    else
        cv_text = [cv_text ' Out-of-sample prediction performance was estimated using the cross-validation classification accuracy for each class, the overall classification accuracy, and the area under the ROC curve. '];   
    end
    if cfg.cv.permutation == 1
        if cfg.cat_Y == 0
            cv_text = [cv_text, ' Permutation testing was performed using ' num2str(cfg.cv.n_perm) ' iterations of the full cross-validation procedure to determine the statistical significance of the cross-validation out-of-sample prediction performance. The average (across folds and repeats) MSE from the cross-validation run was compared against the null distribution obtained from the permutation analyses to determine statistical significance of the full cross-validation results.'];
        else
            cv_text = [cv_text, ' Permutation testing was performed using ' num2str(cfg.cv.n_perm) ' iterations of the full cross-validation procedure to determine the statistical significance of the cross-validation out-of-sample prediction performance. The average (across folds and repeats) area under the ROC curve from the cross-validation run was compared against the null distribution obtained from the permutation analyses to determine statistical significance of the full cross-validtion results.'];
        end
    end  
    if cfg.cat_Y == 1 && isfield(cfg, 'cost')
        if ~isequal(cfg.cost, [0, 1; 1, 0])
            cv_text = [cv_text ' The misclassification cost assigned to the class labeled -1 was set to ' num2str(cfg.cost(2,1)) ', and the misclassification cost assigned to the class labeled 1 was set to ' num2str(cfg.cost(1,2)) '.'];
        end
    end
end

%%% Hyper-parameter optimization methods
if cfg.optimize_hyperparams == 1
    switch cfg.model_spec
        case {'plsr', 'pls_da'}
            param_names = ' k parameter (number of PLS components)';
        case {'lasso', 'ridge', 'logistic_lasso', 'logistic_ridge', 'rlinsvr', 'rlinsvc'}
            param_names = 'regularization parameter (lambda)';
        case {'linsvr', 'kernsvr', 'linsvc', 'kernsvc'}
            param_names = 'box constraint and kernel scale parameter';
        case {'rensemble', 'censemble'}
            if iscellstr(cfg.hp_opt.to_optimize)
                for i = 1:length(cfg.hp_opt.to_optimize)
                    if i == 1
                        param_names = [cfg.hp_opt.to_optimize{i}];
                    elseif i > 1 && i < length(cfg.hp_opt.to_optimize)
                        param_names = [param_names ', ' cfg.hp_opt.to_optimize{i}];
                    else
                        param_names = [param_names ', and ' cfg.hp_opt.to_optimize{i}];
                    end
                end
            elseif isa(cfg.hp_opt.to_optimize(1), 'optimizableVariable')
                for i = 1:length(cfg.hp_opt.to_optimize)
                    if i == 1
                        param_names = [cfg.hp_opt.to_optimize(i).Name];
                    elseif i > 1 && i < length(cfg.hp_opt.to_optimize)
                        param_names = [param_names ', ' cfg.hp_opt.to_optimize(i).Name];
                    else
                        param_names = [param_names ', and ' cfg.hp_opt.to_optimize(i).Name];
                    end
                end
            end
    end
    if cfg.hp_opt.bayes_opt == 0
        if strcmp(cfg.hp_opt.cv_type, 'KFold')
            if cfg.hp_opt.repeats > 1 
                hp_opt_text = [' Hyper-parameter optimization was performed to determine the ' param_names ' using ' num2str(cfg.hp_opt.repeats) ' repetitions of '...
                    num2str(cfg.hp_opt.folds) '-fold cross-validation with repartitioning at each repetition.'];
            else
                hp_opt_text = [' Hyper-parameter optimization was performed to determine the ' param_names ' using ' ...
                    num2str(cfg.hp_opt.folds) '-fold cross-validation with repartitioning at each repetition.'];
            end
        elseif strcmp(cfg.hp_opt.cv_type, 'Holdout')
            if cfg.hp_opt.repeats > 1 
                hp_opt_text = [' Hyper-parameter optimization was performed to determine the ' param_names ' using ' num2str(cfg.hp_opt.repeats) ' repetitions of '...
                    num2str(cfg.hp_opt.holdout) '% holdout cross-validation with repartitioning at each repetition.'];
            else
                hp_opt_text = [' Hyper-parameter optimization was performed to determine the ' param_names ' using ' ...
                    num2str(cfg.hp_opt.holdout .* 100) 'holdout cross-validation with repartitioning at each repetition.'];
            end        
        elseif strcmp(cfg.hp_opt.cv_type, 'LOOCV')
            hp_opt_text = [' Hyper-parameter optimization was performed to determine the ' param_names ' using leave-one-out cross-validation.'];
        end
    else
        if cfg.hp_opt.repartition == 1
            hp_opt_text = [' Hyper-parameter optimization was performed to determine the ' param_names ' via ' cfg.hp_opt.cv_type ' cross-validation using Bayesian optimization with ' num2str(cfg.hp_opt.opt_iter) ' objective iterations and repartitioning at each iteration.'];
        else
            hp_opt_text = [' Hyper-parameter optimization was performed ' param_names ' via ' cfg.hp_opt.cv_type ' cross-validation using Bayesian optimization with ' num2str(cfg.hp_opt.opt_iter) ' objective iterations.'];
        end
    end
end

%%% Inferential modeling methods
if cfg.fit_explanatory_model == 1
    expl_text = [mdl_text, ' An inferential model was fit to the full dataset.'];
    if cfg.optimize_hyperparams == 1 && cfg.cross_validation == 1
        expl_text = [expl_text ' ' hp_opt_text];
    end
    if cfg.permutation == 1
        expl_text = [expl_text ' To determine the significance of the inferential model compared to an empirical null distribution of model fits, permutation testing was performed using ' num2str(cfg.perm.n_perm) ' permutation iterations.'];
        if cfg.perm.coeff_cfwe == 1
            expl_text = [expl_text, ' The continuous FWE method (Mirman et al., 2017) was used to determine a single family-wise error threshold across all predictors at predictor thresholds ' ...
                'of [1,10,50,100,500,1000].'];
        end
        if cfg.perm.coeff_p == 1
            expl_text = [expl_text, ' Coefficient-level permutation p-values were also generated.'];
            if cfg.perm.coeff_fwe == 1
                expl_text = [expl_text, ' Bonferroni-Holm correction was used to obtain FWE-corrected coefficient-level permutation p-values.'];
            end
            if cfg.perm.coeff_fdr == 1
                expl_text = [expl_text, ' The Benjamini-Hochberg procedure was used to obtain FDR-corrected coefficient-level permutation p-values.'];
            end
        end
    end
    if cfg.bootstrap == 1
        if cfg.boot.get_pvals == 1
            expl_text = [expl_text, ' Bootstrap resampling with ' num2str(cfg.boot.n_boot) ' bootstrap iterations was used to estimate the statistical significance of individual coefficients (Kohoutova et al., 2020).'];
            if cfg.boot.coeff_fwe == 1
                expl_text = [expl_text, ' The Bonferroni-Holm procedure was used to obtain FWE-corrected coefficient-level p-values.'];
            end
            if cfg.boot.coeff_fdr == 1
                 expl_text = [expl_text, ' The Benjamini-Hochberg procedure was used to obtain FDR-corrected coefficient-level p-values.'];
            end
        end
        if cfg.boot.get_cis == 1
            if cfg.boot.ci_coeff_fwe == 1 && strcmp(cfg.boot.ci_type, 'normal') 
                expl_text = [expl_text,' ' num2str((1-cfg.boot.ci_alpha_thresh_model).*100) '% bootstrap confidence intervals were estimated for the model fit, and FWE-corrected '  num2str((1-cfg.boot.ci_alpha_thresh_coeff).*100) '% bootstrap confidence intervals were estimated for individual coefficients using the '...
                    'normal approximation method.'];
            else
                expl_text = [expl_text, ' ' num2str((1-cfg.boot.ci_alpha_thresh_model).*100) '% bootstrap confidence intervals were estimated for the model fit, and '  num2str((1-cfg.boot.ci_alpha_thresh_coeff).*100) '% bootstrap confidence intervals were estimated for individual coefficients using the '...
                cfg.boot.ci_type ' method.'];
            end
        end
    end           
    if cfg.cross_validation == 1
        expl_text = [expl_text ' ' cv_text];
    else
        if cfg.cat_Y == 1 && isfield(cfg, 'cost')
            if ~isequal(cfg.cost, [0, 1; 1, 0])
                expl_text = [' The misclassification cost assigned to the class labeled -1 was set to ' num2str(cfg.cost(2,1)) ', and the misclassification cost assigned to the class labeled 1 was set to ' num2str(cfg.cost(1,2)) '.'];
            end
        end
    end
    method_text = expl_text;
elseif cfg.fit_explanatory_model == 0 && cfg.cross_validation == 1
    method_text = [mdl_text, cv_text];
end
