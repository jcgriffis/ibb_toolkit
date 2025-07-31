function [model_results] = fit_and_evaluate_model(cfg)

% Depending on user-defined parameters in the input cfg file, this function can:
% - preprocess input data to remove missing observations and columns with insufficient data
% - generate a lesion frequency map
% - apply direct total lesion volume control to a lesion predictor matrix
% - stratify samples based on a pre-defined stratification variable
% - tune, fit, and evaluate out-of-sample performance for any of the following models through nested cross-validation:
%       - partial least squares regression/classification
%       - ridge, lasso, or linear support vector regression/classification
%       - kernel support vector regression/classification 
% - tune, fit, and evaluate in-sample explanatory power of any of the above models through cross-validation and bootstrapping/permutation tests
% - run mass univariate analyses
% - identify models and predictors with significant predictors through bootstrapping (e.g., see Griffis et al., 2019, Cell Reports; Kohoutova et al., 2020, Nature Protocols).
% - identify significant models and predictors with permutation testing 

% Joseph Griffis, 2023

%%%%%%% Inputs %%%%%%%%
% cfg - structure containing analysis parameters and input data (see documentation)

%%%%%% Outputs %%%%%%%%
% model_results - structure containing relevant results and parameters for fixed effect model and bootstrap analyses
% cv_results - structure containing relevant results for nested CV analyses

% Start timer
tic;

%%%%% Reduce predictor data and check for missing data

% Trim predictor data
disp(['Trimming data to retain columns with at least ' num2str(cfg.freq_thresh) ' observations greater than or equal to ' num2str(cfg.min_obs)]);
[X, X_inc_ind, Y, cfg] = trim_input_data(cfg);

% Remove subjects with no data after trimming and store data check results in results structure
if isfield(cfg, 'cv')
    if ~isfield(cfg.cv, 'stacked_model') || cfg.cv.stacked_model == 0
        [model_results.has_data, X, Y, cfg] = check_for_missing_data(X, Y, cfg);
    else
        model_results.has_data = ones(length(Y),1); % don't drop subjects to allow for model stacking
    end
else
    [model_results.has_data, X, Y, cfg] = check_for_missing_data(X, Y, cfg);
end

% Make sure that Y is appropriate
Y = check_Y_data(X, Y, cfg);

% Get model dimensionality for OLSR
if strcmp(cfg.model_spec, 'olsr')
    mdl = fitlm(X, Y, cfg.model_terms, "CategoricalVars", cfg.categorical_X);
    cfg.mdl_dim = length(mdl.Coefficients.Estimate);
end

% Create output directory if needed
if ~isfolder(cfg.out_dir)
    mkdir(cfg.out_dir);
end
cd(cfg.out_dir);

% Get lesion frequency map
if cfg.gen_freq_map == 1
   get_freq_map(X, X_inc_ind, cfg);
end

% Write out lesion volume correlation map
if cfg.gen_lvol_map == 1
    get_lvol_map(X, X_inc_ind, cfg);
end

% Apply dTLVC to input data if indicated
if cfg.dtlvc == 1
    disp('Applying direct total lesion volume control...');
    X = apply_dtlvc(X);
end

% Stratify if indicated
if ~isempty(cfg.strat_var)
    if cfg.cat_Y == 1 && isequal(cfg.strat_var, Y) % Generate dummy code for categorical outcomes if the grouping variable is the response
        cfg.strat_groups = categorical(cfg.strat_var); % just define stratification groups as a categorical data type version of stratification variable
    elseif ~iscategorical(cfg.strat_var) % if it's already categorical, just leave it as it is, otherwise get quaritle groups 
        cfg.strat_groups = get_quartile_groups(cfg.strat_var); % Get quartile groups for continuous outocmes
    elseif iscategorical(cfg.strat_var) 
        cfg.strat_groups = cfg.strat_var;
    end
end

%%%% Set RNG seed and save
rng(0, 'twister'); % rng seed
rng_seed = {0, 'twister'};

% Run cross-validation if indicated
if cfg.cross_validation == 1
    if ~contains(cfg.model_spec, {'municorr', 'munilr', 'muniolsr', 'bmunz', 'ttest', 'prop_sub', 'munimnr'})
            run_nested_cv(X, Y, cfg, 0);
        if cfg.cv.permutation == 1
            perm_nested_cv(X, Y, cfg);
        end
    else
        warning('Nested CV not supported for mass univariate analysis; skipping to full sample analysis');
    end
end

% Fit explanatory model to full sample (with bootstrapping/permutation testing if indicated)
if cfg.fit_explanatory_model == 1
    if ~isempty(cfg.confounds)
        if isfield(cfg, 'confound_opt')
            if strcmp(cfg.confound_opt, 'regress_both')
                [model_results, X_resid, Y_resid] = fit_explanatory_model(X,Y,model_results,cfg);
                model_results.X_resid = X_resid;
                model_results.Y_resid = Y_resid;
                clear X_resid Y_resid
            end
        elseif ~isfield(cfg, 'confound_opt') && cfg.cat_Y == 0
                [model_results, ~, Y_resid] = fit_explanatory_model(X,Y,model_results,cfg);
                model_results.Y_resid = Y_resid; 
                clear Y_resid
        elseif cfg.cat_Y == 1
                [model_results, X_resid, ~] = fit_explanatory_model(X,Y,model_results,cfg);
                model_results.X_resid = X_resid; 
                clear X_resid
        end
    else
        [model_results, ~, ~] = fit_explanatory_model(X, Y, model_results, cfg);
    end
end

disp('Finished running modeling analysis');
time_to_run = toc;
disp(['Time to to complete:' num2str(round(time_to_run,1) / 60) ' minutes']);

%%%% Store other relevant results in results structure 
model_results.runtime = round(time_to_run,1) / 60;
model_results.X = X; % trimmed predictor matrix
model_results.X_ind = X_inc_ind'; % indices of included cells (i.e. to map back to original predictor matrix)
model_results.rng_seed = rng_seed; % rng seed (for reproducibility)
model_results.cfg = cfg; % append cfg file for reproducibility
if ~isfield(model_results, 'obs_y')
    model_results.obs_y = Y;
end
model_results.methods_text = gen_methods_text(model_results); % summary methods paragraph
model_results.results_text = gen_results_text(model_results); % summary results paragraph

% Save results
cd(cfg.out_dir);

% Write out boilerplate methods and results descriptions
writelines(model_results.methods_text, 'methods_summary.txt');
writelines(model_results.results_text, 'results_summary.txt');

% Save full model results
if cfg.save_model_results == 1
    if isfile('model_results.mat')
       save('model_results.mat', 'model_results', '-append');
    else
       save('model_results.mat', 'model_results', '-v7.3');
    end
end
end