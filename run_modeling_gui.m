classdef run_modeling_gui < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        CatYCheckBox                    matlab.ui.control.CheckBox
        NuisanceRegressorCheckBox       matlab.ui.control.CheckBox
        RegressLesionVolumefromYCheckBox  matlab.ui.control.CheckBox
        SelectNuisanceRegressorsListBox  matlab.ui.control.ListBox
        SelectNuisanceRegressorsListBoxLabel  matlab.ui.control.Label
        SelectOutcomeVariableYListBox   matlab.ui.control.ListBox
        SelectOutcomeVariableYListBoxLabel  matlab.ui.control.Label
        MisclassificationCostsPanel     matlab.ui.container.Panel
        MisClasCostG1                   matlab.ui.control.NumericEditField
        Group1EditField_3Label          matlab.ui.control.Label
        MisClasCostG2                   matlab.ui.control.NumericEditField
        Group1EditField_4Label          matlab.ui.control.Label
        StatisticalThresholdsPanel      matlab.ui.container.Panel
        FDRLabel                        matlab.ui.control.Label
        FDRThresh                       matlab.ui.control.NumericEditField
        UncPThresh                      matlab.ui.control.NumericEditField
        UncEditFieldLabel               matlab.ui.control.Label
        FweThresh                       matlab.ui.control.NumericEditField
        FWEEditFieldLabel               matlab.ui.control.Label
        PermutationTestingPanel         matlab.ui.container.Panel
        RegistryFlag                    matlab.ui.control.CheckBox
        OutDirText                      matlab.ui.control.TextArea
        SelectOutputDirectoryButton     matlab.ui.control.Button
        RunAnalysisButton               matlab.ui.control.Button
        MultivariateModelingToolAnalysisConfigurationLabel  matlab.ui.control.Label
        UseParallelProcessingRecommendedCheckBox  matlab.ui.control.CheckBox
        FitExplanatoryModelCheckBox     matlab.ui.control.CheckBox
        ResultOptionsPanel              matlab.ui.container.Panel
        WritePermutationImagesCheckBox  matlab.ui.control.CheckBox
        WriteBootstrapImagesCheckBox    matlab.ui.control.CheckBox
        SelectModelDropDown             matlab.ui.control.DropDown
        SelectModelLabel                matlab.ui.control.Label
        BrainMaskText                   matlab.ui.control.TextArea
        ImDirText                       matlab.ui.control.TextArea
        CSVTtext                        matlab.ui.control.TextArea
        PredictorTypeButtonGroup        matlab.ui.container.ButtonGroup
        OtherPredButton                 matlab.ui.control.RadioButton
        LesionButton                    matlab.ui.control.RadioButton
        MinimumDataThresholdsPanel      matlab.ui.container.Panel
        MinObs                          matlab.ui.control.NumericEditField
        MinimumObservationLabel         matlab.ui.control.Label
        MinFreq                         matlab.ui.control.NumericEditField
        MinimumFrequencyLabel           matlab.ui.control.Label
        SelectCSVfileButton             matlab.ui.control.Button
        SelectDataDirectoryButton       matlab.ui.control.Button
        HyperParameterTuningPanel       matlab.ui.container.Panel
        OptimizeHyperparametersCheckBox  matlab.ui.control.CheckBox
        HpOptRepeats                    matlab.ui.control.NumericEditField
        HpOptFolds                      matlab.ui.control.NumericEditField
        HpOptCvType                     matlab.ui.control.DropDown
        HpOptNumberofRepeatsLabel       matlab.ui.control.Label
        HpOptNumberofFoldsLabel         matlab.ui.control.Label
        HpOptCVStrategyLabel            matlab.ui.control.Label
        SelectBrainMaskButton           matlab.ui.control.Button
        PredictorModalityButtonGroup    matlab.ui.container.ButtonGroup
        ArrayButton                     matlab.ui.control.RadioButton
        MatrixButton                    matlab.ui.control.RadioButton
        NIFTIButton                     matlab.ui.control.RadioButton
        CrossValidationOptionsNestedifHyperParameterTuningEnabledPanel  matlab.ui.container.Panel
        EnableModelStackingCheckBox     matlab.ui.control.CheckBox
        CVPermIter                      matlab.ui.control.NumericEditField
        PermutationsEditField_2Label    matlab.ui.control.Label
        CVPermPanel                     matlab.ui.container.Panel
        CVPermCheckBox                  matlab.ui.control.CheckBox
        RunCrossValidationCheckBox      matlab.ui.control.CheckBox
        CrossValidationSettingsPanel    matlab.ui.container.Panel
        OuterRepeats                    matlab.ui.control.NumericEditField
        OuterFolds                      matlab.ui.control.NumericEditField
        CVType                          matlab.ui.control.DropDown
        NumberofRepeatsLabel_2          matlab.ui.control.Label
        NumberofOuterFoldsLabel         matlab.ui.control.Label
        CVStrategyLabel_2               matlab.ui.control.Label
        GeneralModelingOptionsPanel     matlab.ui.container.Panel
        BootIter                        matlab.ui.control.NumericEditField
        BootAnalysesPanel               matlab.ui.container.Panel
        BootstrapLabel                  matlab.ui.control.Label
        PermutationsEditFieldLabel      matlab.ui.control.Label
        BootCheckBox                    matlab.ui.control.CheckBox
        PermCheckBox                    matlab.ui.control.CheckBox
        StratifyTrainTestSplitsonYCheckBox  matlab.ui.control.CheckBox
        ConfIntWidth                    matlab.ui.control.NumericEditField
        PermIter                        matlab.ui.control.NumericEditField
        ConfIntLabel                    matlab.ui.control.Label
        BootCIsCheckBox                 matlab.ui.control.CheckBox
        PermVoxPValsCheckBox            matlab.ui.control.CheckBox
        PermApplycFWECheckBox           matlab.ui.control.CheckBox
        ApplyDTLVCCheckBox              matlab.ui.control.CheckBox
        StandardizeXCheckBox            matlab.ui.control.CheckBox
        StandardizeYCheckBox            matlab.ui.control.CheckBox
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            model_dir = fileparts(which('run_modeling_gui.mlapp'));
            addpath(model_dir);
        end

        % Callback function: CSVTtext, SelectCSVfileButton
        function SelectCSVfileButtonPushed(app, event)
            model_dir = fileparts(which('run_modeling_gui.mlapp'));            
            if isfile(fullfile(model_dir, 'default_paths.mat'))
                load(fullfile(model_dir, 'default_paths.mat'), 'csv_dir');
                [csv_file, csv_path] = uigetfile({'*.csv'; '*xls'; '*xlsx'; '*.'},...
                    'Select CSV file containing behavioral data', csv_dir);
            else
                [csv_file, csv_path] = uigetfile({'*.csv'; '*xls'; '*xlsx'; '*.'},...
                    'Select CSV file containing behavioral data');  
            end
            if csv_file == 0
                error('No file selected; please try again');
            else
                app.CSVTtext.Value = fullfile(csv_path, csv_file);
            end
            my_csv = readtable(fullfile(csv_path, csv_file));
            my_vars = my_csv.Properties.VariableNames;
            my_vars = setdiff(my_vars, "study_id");
            app.SelectOutcomeVariableYListBox.Items = my_vars;
            if app.NuisanceRegressorCheckBox.Value == 1
                app.SelectNuisanceRegressorsListBox.Items = my_vars;
            end
        end

        % Button pushed function: SelectDataDirectoryButton
        function SelectDataDirectoryButtonPushed(app, event)
            model_dir = fileparts(which('run_modeling_gui.mlapp'));            
            if isfile(fullfile(model_dir, 'default_paths.mat'))
                load(fullfile(model_dir, 'default_paths.mat'), 'data_dir');
                cfg.img_dir = uigetdir(data_dir, 'Select directory containing imaging data');
            else
                cfg.img_dir = uigetdir(pwd, 'Select directory containing imaging data');
            end
            if cfg.img_dir == 0
                error('No directory selected; please try again');
            else
                app.ImDirText.Value = cfg.img_dir;
            end           
        end

        % Button pushed function: SelectBrainMaskButton
        function SelectBrainMaskButtonPushed(app, event)
            model_dir = fileparts(which('run_modeling_gui.mlapp'));            
            if app.NIFTIButton.Value == 1
                if isfile(fullfile(model_dir, 'default_paths.mat'))
                    load(fullfile(model_dir, 'default_paths.mat'), 'mask_file');
                    [mask_file, mask_path] = uigetfile({'*.nii.gz'; '*.nii'; '*.'}, 'Select brain mask file',...
                        mask_file);
                else
                    [mask_file, mask_path] = uigetfile({'*.nii.gz'; '*.nii'; '*.'}, 'Select brain mask file', pwd);
                end
                if mask_file == 0
                    error('No file selected; please try again');
                else
                    app.BrainMaskText.Value = fullfile(mask_path, mask_file);
                end
            elseif app.MatrixButton.Value == 1
                [parcel_file, parcel_path] = uigetfile({'*.mat'}, 'Select parcellation_table',...
                    pwd);       
                if parcel_file == 0
                    error('No file selected; please try again');
                else
                    app.BrainMaskText.Value = fullfile(parcel_path, parcel_file);
                end            
            end
        end

        % Button pushed function: RunAnalysisButton
        function RunAnalysisButtonPushed(app, event)
            
            % Get model directory
            model_dir = fileparts(which('run_modeling_gui.mlapp'));
            cfg.model_dir = fullfile(model_dir, 'model_code');
            
            % Add core functions to path
            addpath(fullfile(cfg.model_dir));
            addpath(genpath(fullfile(cfg.model_dir, 'common')));
            
            % Set analysis modality
            if app.LesionButton.Value == 1
                cfg.modality = 'lesion';
            else
                cfg.modality = 'other';
            end
            
            % Set imaging flag
            if app.NIFTIButton.Value == 1
                cfg.write_nifti_images = 1;
                cfg.mask_path = string(app.BrainMaskText.Value);
            elseif app.TableButton.Value == 1
                cfg.write_nifti_images = 0;
            elseif app.MatrixButton.Value == 1
                cfg.write_nifti_images = 0;
                cfg.parcel_table = string(app.BrainMaskText.Value);
            elseif app.ArrayButton.Value == 1
                cfg.write_nifti_images = 0;
            end

            % Set up basic cfg file
            cfg.model_spec = string(app.SelectModelDropDown.Value);
            cfg.fit_explanatory_model = app.FitExplanatoryModelCheckBox.Value;
            cfg.bootstrap = app.BootCheckBox.Value;
            cfg.jackknife = 0;
            cfg.permutation = app.PermCheckBox.Value;
            cfg.optimize_hyperparams = app.OptimizeHyperparametersCheckBox.Value;
            cfg.cross_validation = app.RunCrossValidationCheckBox.Value;         
            cfg.out_dir = string(app.OutDirText.Value);
            cfg.stacked_model = app.EnableModelStackingCheckBox.Value;

            % Get default cfg
            cfg = get_default_model_cfg(cfg);

            %%%%%% Update defaults with user choices

            %%%% General analysis params

            % Data thresholds
            cfg.min_obs = app.MinObs.Value;
            cfg.freq_thresh = app.MinFreq.Value;

            % Standardization
            if app.StandardizeXCheckBox.Value == 1 && app.StandardizeYCheckBox.Value == 1
                cfg.standardize = 3;
            elseif app.StandardizeYCheckBox.Value == 1 && app.StandardizeXCheckBox.Value == 0
                cfg.standardize = 2;
            elseif app.StandardizeXCheckBox.Value == 1 && app.StandardizeYCheckBox.Value == 0
                cfg.standardize = 1;
            else
                cfg.standardize = 0;
            end
            
            % DTLVC 
            cfg.dtlvc = app.ApplyDTLVCCheckBox.Value;

            % Hyperparameter optimization
            cfg.hp_opt.cv_type = string(app.HpOptCvType.Value);        
            cfg.hp_opt.folds = app.HpOptFolds.Value;
            cfg.hp_opt.repeats = app.HpOptRepeats.Value;           
            
            % Statistical testing
            cfg.boot.n_boot = app.BootIter.Value;
            cfg.perm.n_perm = app.PermIter.Value;
            if app.UncPThresh.Value ~= 0
                cfg.unc_thresh = app.UncPThresh.Value;
            else
                cfg.unc_thresh = 0.001;
            end
            if app.FDRThresh.Value ~= 0
                cfg.fdr_thresh = app.FDRThresh.Value;
            else
                cfg.fdr_thresh = 0.05;
            end
            if app.FweThresh.Value ~= 0
                cfg.fwe_thresh = app.FweThresh.Value;
            else
                cfg.fwe_thresh = 0.05;
            end
            cfg.perm.coeff_p = app.PermVoxPValsCheckBox.Value;
            if cfg.perm.coeff_p == 1
                cfg.perm.coeff_fdr = 1;
                cfg.perm.coeff_fwe = 1;
            end
            cfg.perm.coeff_cfwe = app.PermApplycFWECheckBox.Value;

            % Cross-validation
            if app.EnableModelStackingCheckBox == 1
                cfg.cv.stacked_model = 1;
            else
                cfg.cv.stacked_model = 0;
            end
            cfg.cv.cv_type = char(app.CVType.Value);
            cfg.cv.folds = app.OuterFolds.Value;
            cfg.cv.repeats = app.OuterRepeats.Value;
            cfg.cv.permutation = app.CVPermCheckBox.Value;
            cfg.cv.n_perm = app.CVPermIter.Value;

            % Parallel processing
            cfg.parallel = app.UseParallelProcessingRecommendedCheckBox.Value;
            
            % Image output
            cfg.write_perm_images = app.WritePermutationImagesCheckBox.Value;
            cfg.write_boot_images = app.WriteBootstrapImagesCheckBox.Value;
            cfg.write_jack_images = 0;

            % Format and load behavioral and image data 
            beh_csv = readtable(string(app.CSVTtext.Value));
            beh_col = app.SelectOutcomeVariableYListBox.Value;
            id_col = "study_id";
            
            % Get data and put into cfg as long as predictors are NIFTI lesions
            if app.NIFTIButton.Value == 1 
                data_dir = app.ImDirText.Value;
                if app.LesionButton.Value == 1
                    cfg = get_and_format_lesion_data(beh_csv, id_col, beh_col, app.RegistryFlag.Value, data_dir, cfg);
                else
                    cfg = get_and_format_nifti_data(beh_csv, id_col, beh_col, app.RegistryFlag.Value, data_dir, cfg);
                end
            elseif app.MatrixButton.Value == 1
                data_dir = app.ImDirText.Value;
                cfg = get_and_format_matrix_data(beh_csv, id_col, beh_col, app.RegistryFlag.Value, data_dir, cfg);
            elseif app.ArrayButton.Value == 1
                data_dir = app.ImDirText.Value;
                cfg = get_and_format_array_data(beh_csv, id_col, beh_col, app.RegistryFlag.Value, data_dir, cfg);              
            end

            % Set up confounds for nuisance regression if indicated
            disp('Nuisance Regressors: ');
            cfg.confounds = [];
            for i = 1:length(app.SelectNuisanceRegressorsListBox.Value)
                disp(app.SelectNuisanceRegressorsListBox.Value{i});
                cfg.confounds(:,i) = beh_csv.(app.SelectNuisanceRegressorsListBox.Value{i});
                cfg.confound_names{i} = app.SelectNuisanceRegressorsListBox.Value{i};
            end
            if ~isempty(cfg.confounds)
                cfg.confounds = cfg.confounds(cfg.include_inds,:);
                if ~isempty(find(isnan(cfg.confounds)))
                    disp('Checking nuisance regressors...')
                    disp('Some subjects have missing data for nuisance regressors; removing subjects with missing data...');                
                    [cfg.confounds, cfg.missing_inds] = rmmissing(cfg.confounds, 1);
                    cfg.X = cfg.X(cfg.missing_inds == 0, :);
                    cfg.Y = cfg.Y(cfg.missing_inds == 0, :);
                    if isfield(cfg, 'lvol')
                        cfg.lvol = cfg.lvol(cfg.missing_inds == 0);
                    end
                    disp([num2str(numel(find(cfg.missing_inds==1))) ' subjects were removed due to missing data for nuisance regressors.']);
                    disp('Removed subjects will have values of 1 in model_results.cfg.missing_inds');
                end
            end
                
            if app.LesionButton.Value == 1 && app.RegressLesionVolumefromYCheckBox.Value == 1
                cfg.confounds = [cfg.confounds, zscore(cfg.lvol)];
            else
                cfg.confounds = cfg.confounds;
            end

            % Do stratification on outcome if indicated
            if app.StratifyTrainTestSplitsonYCheckBox.Value == 1
                cfg.strat_var = cfg.Y;
            else
                cfg.strat_var = [];
            end
            
            % Set misclassification costs if relevant
            if cfg.cat_Y==1
                cfg.cost(1,2) = app.MisClasCostG1.Value;
                cfg.cost(2,1) = app.MisClasCostG2.Value;
            end
            
            % Set grouping flag for T-test
            if strcmp(cfg.model_spec, 'ttest')
                if app.CatYCheckBox.Value == 1
                    cfg.group_var = 'Y';
                    cfg.cat_Y = 1;
                else
                    cfg.group_var = 'X';
                    cfg.cat_Y = 0;
                end
            end

            % Train and evaluate model
            model_results = fit_and_evaluate_model(cfg);

            % Write out images if indicated
            if app.NIFTIButton.Value == 1
                write_model_weight_maps(model_results);
            elseif app.MatrixButton.Value == 1
                write_model_weights_connectome(model_results);
            end
            
            % Print main results to command window and save to text file
            cd(cfg.out_dir);

            % Print results for regression models 
            if ~contains(cfg.model_spec, ["municorr", "bmunz", "ttest", "munilr", "muniolsr", "prop_sub"])
                if cfg.fit_explanatory_model == 1 && cfg.cat_Y == 0
                    
                   % Print model R-squared
                   diary('explanatory_model_results.txt')
                   disp('--------Model Result Summary-------');
                   disp('-----Inferential Model Results-----');
                   disp('-----------------------------------');
                   disp(['Model R-squared: ' num2str(model_results.r2)]);
                    
                   % Print additional statistics depending on analysis settings
                   if cfg.bootstrap == 1 && cfg.boot.get_cis == 1
                       disp([num2str(100-(10*cfg.boot.ci_alpha_thresh_model)) '% Bootstrap CIs on R-squared:']);
                       disp(['Lower: ' num2str(model_results.boot.r2_ci(1))]);
                       disp(['Upper: ' num2str(model_results.boot.r2_ci(2))]);
                   end
                   if cfg.permutation == 1 
                       disp(['Permutation Test: p=' num2str(model_results.perm.model_pval)]);
                   end
                   if cfg.cross_validation == 1
                       load(fullfile(cfg.out_dir, 'cv_results.mat'));
                       disp(['Cross-validation Correlation: r=' num2str(cv_results.all.corr) ', p=' num2str(cv_results.all.corr_pval)]);
                       if isfield(cv_results, 'perm_pval')
                           disp('-----Predictive Model Results-----');
                           disp(['Cross-Validation Permutation Test: p= ' num2str(cv_results.perm_pval)]);
                           disp(['Mean cross-validation R-squared: ' num2str(mean(cv_results.avg.r2_ss(:)))]);
                       end
                   end
                   diary('off');
                elseif cfg.fit_explanatory_model == 0 && cfg.cat_Y == 0 && cfg.cross_validation == 1
                       load(fullfile(cfg.out_dir, 'cv_results.mat'));   
                       if isfield(cv_results, 'perm_pval')
                           disp('-----Predictive Model Results-----');                       
                           disp(['Cross-Validation Permutation Test: p= ' num2str(cv_results.perm_pval)]);
                           disp(['Mean cross-validation R-squared: ' num2str(mean(cv_results.avg.r2_ss(:)))]);
                       else
                           disp('-----Predictive Model Results-----');                       
                           disp(['Mean cross-validation R-squared: ' num2str(mean(cv_results.avg.r2_ss(:)))]);                      
                       end
                end            
                
                % Print results for classification models 
                if cfg.fit_explanatory_model == 1 && cfg.cat_Y == 1
                    
                   % Print model R-squared
                   diary('explanatory_model_results.txt')
                   disp('--------Model Result Summary-------');
                   disp('-----Inferential Model Results-----');
                   disp('-----------------------------------')
                   disp(['Model Accuracy (Group -1): ' num2str(model_results.classrate(1))]);
                   disp(['Model Accuracy (Group 1): ' num2str(model_results.classrate(2))]);               
                   disp(['Model Accuracy (Overall): ' num2str(model_results.classrate(3))]);
                   disp(['Model Area Under ROC Curve: ' num2str(model_results.roc_auc)]);
     
                   % Print additional statistics depending on analysis settings
                   if cfg.bootstrap == 1 && cfg.boot.get_cis == 1
                       disp([num2str(100-(10*cfg.boot.ci_alpha_thresh_model)) '% Bootstrap CIs on Group -1 Classification Rate:']);
                       disp(['Lower: ' num2str(model_results.boot.classrate_ci(1,1))]);
                       disp(['Upper: ' num2str(model_results.boot.classrate_ci(2,1))]);
                       disp([num2str(100-(10*cfg.boot.ci_alpha_thresh_model)) '% Bootstrap CIs on Group 1 Classification Rate:']);
                       disp(['Lower: ' num2str(model_results.boot.classrate_ci(1,2))]);
                       disp(['Upper: ' num2str(model_results.boot.classrate_ci(2,2))]);
                       disp([num2str(100-(10*cfg.boot.ci_alpha_thresh_model)) '% Bootstrap CIs on Overall Classification Rate:']);
                       disp(['Lower: ' num2str(model_results.boot.classrate_ci(1,3))]);
                       disp(['Upper: ' num2str(model_results.boot.classrate_ci(2,3))]);
                       disp([num2str(100-(10*cfg.boot.ci_alpha_thresh_model)) '% Bootstrap CIs on ROC Area Under Curve:']);
                       disp(['Lower: ' num2str(model_results.boot.roc_auc_ci(1,1))]);
                       disp(['Upper: ' num2str(model_results.boot.roc_auc_ci(2,1))]);                   
                   end
                   if cfg.permutation == 1 
                      disp(['Permutation Test (AUC): p=' num2str(model_results.perm.model_pval)]);                   
                   end
                   if cfg.cross_validation == 1
                      load(fullfile(cfg.out_dir, 'cv_results.mat'));
                      disp(['Cross-validation Fisher Exact Test: stat=' num2str(cv_results.all.fisher_stat.OddsRatio) ', p=' num2str(cv_results.all.fisher_pval)]);
                      if isfield(cv_results, 'perm_pval')
                          disp('-----Predictive Model Results-----');
                          disp(['Cross-Validation Permutation Test: p= ' num2str(cv_results.perm_pval)]);
                          disp(['Mean Cross-Validation Area Under ROC Curve: ' num2str(mean(cv_results.avg.roc_auc(:)))]);             
                      else
                          disp('-----Predictive Model Results-----');
                          disp(['Mean Cross-Validation Area Under ROC Curve: ' num2str(mean(cv_results.avg.roc_auc(:)))]);    
                      end
                   end
                   diary('off');
                elseif cfg.fit_explanatory_model == 0 && cfg.cat_Y == 1 && cfg.cross_validation == 1
                      load(fullfile(cfg.out_dir, 'cv_results.mat'));
                      if isfield(cv_results, 'perm_pval')
                          disp('-----Predictive Model Results-----');
                          disp(['Cross-Validation Permutation Test: p= ' num2str(cv_results.perm_pval)]);
                          disp(['Mean Cross-Validation Area Under ROC Curve: ' num2str(mean(cv_results.avg.roc_auc(:)))]);             
                      else
                          disp('-----Predictive Model Results-----');
                          disp(['Mean Cross-Validation Area Under ROC Curve: ' num2str(mean(cv_results.avg.roc_auc(:)))]);    
                      end
                end
            end

        end

        % Button pushed function: SelectOutputDirectoryButton
        function SelectOutputDirectoryButtonPushed(app, event)
             [out_dir] = uigetdir(pwd, 'Select directory to save results');
            if out_dir == 0
                error('No directory selected; please try again');
            else
                app.OutDirText.Value = out_dir;
            end           
        end

        % Value changed function: HpOptCvType
        function HpOptCvTypeValueChanged(app, event)
            value = app.HpOptCvType.Value;
            if strcmp(string(value), 'LOOCV')
                app.HpOptFolds.Enable = 'off';
                app.HpOptFolds.Value = 'N/A';
                app.HpOptRepeats.Value = 1;
            else
                app.HpOptFolds.Enable = 'on';
                app.HpOptFolds.Value = 5;
                app.HpOptRepeats.Value = 5;
            end
        end

        % Value changed function: StratifyTrainTestSplitsonYCheckBox
        function StratifyTrainTestSplitsonYCheckBoxValueChanged(app, event)

           
        end

        % Value changed function: PermCheckBox
        function PermCheckBoxValueChanged(app, event)
            value = app.PermCheckBox.Value;
            if value == 0
                app.PermIter.Enable = 'off';
                app.PermIter.Value = 0;
                app.PermApplycFWECheckBox.Enable = 'off';
                app.PermApplycFWECheckBox.Value = 0;
                if app.BootCheckBox.Value == 0                
                    app.FweThresh.Value = 0.05;
                    app.FweThresh.Enable = 'off';
                    app.FDRThresh.Value = 0.05;
                    app.FDRThresh.Enable = 'off';
                    app.UncPThresh.Value = 0.001;
                    app.UncPThresh.Enable = 'off';       
                else
                    app.FweThresh.Value = 0.05;
                    app.FweThresh.Enable = 'on';
                    app.FDRThresh.Value = 0.05;
                    app.FDRThresh.Enable = 'on';
                    app.UncPThresh.Value = 0.001;
                    app.UncPThresh.Enable = 'on';
                end                
            else
                app.PermIter.Enable = 'on';
                app.PermIter.Value = 1000;
                app.PermApplycFWECheckBox.Enable = 'on';
                app.PermApplycFWECheckBox.Value = 1;    
                app.FweThresh.Value = 0.05;
                app.FweThresh.Enable = 'on';
                app.FDRThresh.Value = 0.05;
                app.FDRThresh.Enable = 'on';
                app.UncPThresh.Value = 0.001;
                app.UncPThresh.Enable = 'on';                
            end
        end

        % Value changed function: BootCheckBox
        function BootCheckBoxValueChanged(app, event)
            value = app.BootCheckBox.Value;
            if value == 0
                app.BootIter.Enable = 'off';
                app.BootIter.Value = 0;
                app.BootCIsCheckBox.Enable = "off";
                app.BootCIsCheckBox.Value = 0;
                app.ConfIntWidth.Value = 0;
                app.ConfIntWidth.Enable = 'off';
                if app.PermCheckBox.Value == 0                
                    app.FweThresh.Value = 0.05;
                    app.FweThresh.Enable = 'off';
                    app.FDRThresh.Value = 0.05;
                    app.FDRThresh.Enable = 'off';
                    app.UncPThresh.Value = 0.001;
                    app.UncPThresh.Enable = 'off';       
                else
                    app.FweThresh.Value = 0.05;
                    app.FweThresh.Enable = 'on';
                    app.FDRThresh.Value = 0.05;
                    app.FDRThresh.Enable = 'on';
                    app.UncPThresh.Value = 0.001;
                    app.UncPThresh.Enable = 'on';
                end
            else
                app.BootIter.Enable = 'on';
                app.BootIter.Value = 1000;
                app.BootCIsCheckBox.Enable = 'on';
                app.BootCIsCheckBox.Value = 1;   
                app.ConfIntWidth.Value = 95;
                app.ConfIntWidth.Enable = 'on';    
                app.FweThresh.Value = 0.05;
                app.FweThresh.Enable = 'on';
                app.FDRThresh.Value = 0.05;
                app.FDRThresh.Enable = 'on';
                app.UncPThresh.Value = 0.001;
                app.UncPThresh.Enable = 'on';
            end            
        end

        % Selection changed function: PredictorTypeButtonGroup
        function PredictorTypeButtonGroupSelectionChanged(app, event)
            if app.LesionButton.Value == 1 && ~contains(app.SelectModelDropDown.Value, ["ttest", "bmunz", "municorr", "munilr"])
                app.ApplyDTLVCCheckBox.Enable = "on";
                app.ApplyDTLVCCheckBox.Value = 1;
                app.RegressLesionVolumefromYCheckBox.Enable = "on";
                app.RegressLesionVolumefromYCheckBox.Value = 0;
                app.MinFreq.Value = 4;
                app.MinObs.Value = 1;                
            elseif app.LesionButton.Value == 1 && contains(app.SelectModelDropDown.Value, ["ttest", "bmunz", "municorr", "munilr"])
                app.ApplyDTLVCCheckBox.Enable = "on";
                app.ApplyDTLVCCheckBox.Value = 0;
                app.RegressLesionVolumefromYCheckBox.Enable = "on";
                app.RegressLesionVolumefromYCheckBox.Value = 1;
                app.MinFreq.Value = 4;
                app.MinObs.Value = 1;                     
            else
                app.ApplyDTLVCCheckBox.Value = 0;
                app.ApplyDTLVCCheckBox.Enable = 'off';
                app.RegressLesionVolumefromYCheckBox.Enable = "off";
                app.RegressLesionVolumefromYCheckBox.Value = 0;
                app.MinFreq.Value = 1;
                app.MinObs.Value = 0;
            end
        end

        % Value changed function: RunCrossValidationCheckBox
        function RunCrossValidationCheckBoxValueChanged(app, event)
            value = app.RunCrossValidationCheckBox.Value;
            if value == 1
                app.CVType.Enable = "on";
                app.CVType.Value = "KFold";
                app.CVPermIter.Enable = "on";
                app.CVPermIter.Value = 20;
                app.OuterFolds.Enable = "on";
                app.OuterFolds.Value = 5;
                app.OuterRepeats.Enable = "on";
                app.OuterRepeats.Value = 5;
                app.CVPermCheckBox.Enable = "on";
                app.CVPermIter.Enable = "on";
            else
                app.CVType.Enable = "off";
                app.CVType.Value = "KFold";
                app.CVPermIter.Enable = "off";
                app.CVPermIter.Value = 0;
                app.OuterFolds.Enable = "off";
                app.OuterFolds.Value = 0;
                app.OuterRepeats.Enable = "off";
                app.OuterRepeats.Value = 0;
                app.CVPermCheckBox.Enable = "off";
                app.CVPermIter.Enable = "off";
            end
        end

        % Value changed function: SelectModelDropDown
        function SelectModelDropDownValueChanged(app, event)
            value = app.SelectModelDropDown.Value;
            if ~contains(value, ["logistic", "pls_da", "svc", "censemble"])
                app.MisClasCostG1.Enable = 0;
                app.MisClasCostG2.Enable = 0;
            else
                app.MisClasCostG1.Enable = 1;
                app.MisClasCostG2.Enable = 1;
            end
            if contains(value, ["logistic", "pls_da", "svc", "censemble", "munilr"])
                app.RegressLesionVolumefromYCheckBox.Text = 'Regress Lesion Volume from X';
            else
                app.RegressLesionVolumefromYCheckBox.Text = 'Regress Lesion Volume from Y';
            end
            if contains(value, ["municorr", "ttest", "bmunz", "munilr", "muniolsr", ...
                    'prop_sub'])
                % Bootstrap options
                app.BootCheckBox.Enable = 0;
                app.BootCIsCheckBox.Enable = 0;
                app.BootstrapLabel.Enable = 0;
                app.BootCIsCheckBox.Enable = 0;
                app.BootIter.Enable = 0;
                app.BootIter.Value = 0;
                app.BootCheckBox.Value = 0;
                app.BootCIsCheckBox.Value = 0;   
                app.ConfIntWidth.Enable = 0;                
                app.WriteBootstrapImagesCheckBox.Enable = 0;
                app.WriteBootstrapImagesCheckBox.Value = 0;
                % Hyperparameter optimization options
                app.HpOptRepeats.Value = 0;
                app.HpOptFolds.Value = 0;
                app.OptimizeHyperparametersCheckBox.Value = 0;
                app.OptimizeHyperparametersCheckBox.Enable = 0;
                app.HpOptCvType.Enable = 0;
                app.HpOptFolds.Enable = 0;
                app.HpOptRepeats.Enable = 0;
                % Permutation testing options
                app.PermCheckBox.Value = 1;
                app.PermCheckBox.Enable = 1;                    
                app.PermApplycFWECheckBox.Value = 1;
                app.PermApplycFWECheckBox.Enable = 1;                    
                app.PermIter.Value = 10000;
                app.PermIter.Enable = 1;           
                app.FweThresh.Value = 0.05;
                app.FweThresh.Enable = 1;
                app.WritePermutationImagesCheckBox.Value = 1;
                app.WritePermutationImagesCheckBox.Enable = 1;                
                if ~strcmp(value, 'prop_sub')
                    app.PermVoxPValsCheckBox.Value = 0;
                    app.PermVoxPValsCheckBox.Enable = 1;
                    app.StandardizeXCheckBox.Enable = 1;
                    app.StandardizeXCheckBox.Value = 0;
                    app.StandardizeYCheckBox.Enable = 1;
                    app.StandardizeYCheckBox.Value = 0;
                    app.FDRThresh.Value = 0.05;
                    app.FDRThresh.Enable = 1;
                    app.UncPThresh.Value = .001;
                    app.UncPThresh.Enable = 1;      
                    app.ApplyDTLVCCheckBox.Enable = 1;
                    app.ApplyDTLVCCheckBox.Value = 0;
                    app.MinFreq.Value = 10;
                    app.RegressLesionVolumefromYCheckBox.Enable = 1;
                    app.RegressLesionVolumefromYCheckBox.Value = 1;
                    app.SelectNuisanceRegressorsListBox.Enable = 1;
                    app.NuisanceRegressorCheckBox.Enable = 1;                    
                else
                    app.PermVoxPValsCheckBox.Value = 0;
                    app.PermVoxPValsCheckBox.Enable = 0;
                    app.StandardizeXCheckBox.Enable = 0;
                    app.StandardizeXCheckBox.Value = 0;
                    app.StandardizeYCheckBox.Enable = 0;
                    app.StandardizeYCheckBox.Value = 0;
                    app.FDRThresh.Value = 0;
                    app.FDRThresh.Enable = 0;
                    app.FweThresh.Value = 0;
                    app.FweThresh.Enable = 0;
                    app.UncPThresh.Value = 0;
                    app.UncPThresh.Enable = 0;
                    app.ApplyDTLVCCheckBox.Enable = 0;
                    app.ApplyDTLVCCheckBox.Value = 0;     
                    app.MinFreq.Value = 1;
                    app.RegressLesionVolumefromYCheckBox.Enable = 0;
                    app.RegressLesionVolumefromYCheckBox.Value = 0;
                    app.SelectNuisanceRegressorsListBox.Enable = 0;
                    app.SelectNuisanceRegressorsListBox.Value = {};
                    app.NuisanceRegressorCheckBox.Value = 0;
                    app.NuisanceRegressorCheckBox.Enable = 0;
                end
                % Cross-validation options
                app.CVType.Enable = 0;
                app.CVPermPanel.Enable = 0;
                app.EnableModelStackingCheckBox.Enable = 0;
                app.StratifyTrainTestSplitsonYCheckBox.Enable = 0;
                app.RunCrossValidationCheckBox.Enable = 0;    
                app.RunCrossValidationCheckBox.Value = 0;
                app.CVPermIter.Value = 0;
                app.CVPermIter.Enable = 0;
                app.EnableModelStackingCheckBox.Value = 0;
                app.OuterFolds.Value = 0;
                app.OuterFolds.Enable = 0;
                app.OuterRepeats.Value = 0;
                app.OuterRepeats.Enable = 0;
                % General options
                app.ApplyDTLVCCheckBox.Value = 0;
                % T-test options
                if strcmp(value, 'ttest')
                    app.CatYCheckBox.Visible = 1;
                    app.CatYCheckBox.Enable = 1;
                    app.CatYCheckBox.Value = 1;
                else
                    app.CatYCheckBox.Visible = 0;
                    app.CatYCheckBox.Enable = 0;
                    app.CatYCheckBox.Value = 0;
                end
            else
                % Bootstrap options
                app.BootCheckBox.Enable = 1;
                app.BootCIsCheckBox.Enable = 1;
                app.BootstrapLabel.Enable = 1;
                app.BootCIsCheckBox.Enable = 1;
                app.BootIter.Enable = 1;
                app.BootIter.Value = 1000;
                app.BootCheckBox.Value = 1;
                app.BootCIsCheckBox.Value = 1;   
                app.ConfIntWidth.Enable = 1;                
                app.WriteBootstrapImagesCheckBox.Enable = 1;
                app.WriteBootstrapImagesCheckBox.Value = 1;
                % Hyperparameter optimization options
                app.HpOptRepeats.Value = 5;
                app.HpOptFolds.Value = 5;
                app.OptimizeHyperparametersCheckBox.Value = 1;
                app.OptimizeHyperparametersCheckBox.Enable = 1;
                app.HpOptCvType.Enable = 1;
                app.HpOptFolds.Enable = 1;
                app.HpOptRepeats.Enable = 1;
                % Permutation testing options
                app.PermCheckBox.Value = 0;
                app.PermCheckBox.Enable = 1;                                    
                app.PermApplycFWECheckBox.Value = 0;
                app.PermApplycFWECheckBox.Enable = 1;                                    
                app.PermIter.Value = 0;
                app.PermIter.Enable = 1;
                app.PermVoxPValsCheckBox.Value = 0;
                app.PermVoxPValsCheckBox.Enable = 1;
                app.WritePermutationImagesCheckBox.Value = 0;
                app.WritePermutationImagesCheckBox.Enable = 1;             
                % Cross-validation options
                app.CVType.Enable = 1;
                app.CVPermPanel.Enable = 1;
                app.EnableModelStackingCheckBox.Enable = 1;
                app.StratifyTrainTestSplitsonYCheckBox.Enable = 1;
                app.RunCrossValidationCheckBox.Enable = 1;    
                app.RunCrossValidationCheckBox.Value = 1;
                app.CVPermIter.Value = 50;
                app.CVPermIter.Enable = 1;
                app.OuterFolds.Enable = 1;
                app.OuterFolds.Value = 5;
                app.OuterRepeats.Enable = 1;
                app.OuterRepeats.Value = 5;
                app.EnableModelStackingCheckBox.Value = 1;
                % T-test options
                app.CatYCheckBox.Visible = 0;
                app.CatYCheckBox.Enable = 0;
                app.CatYCheckBox.Value = 0;
                % Other
                app.StandardizeXCheckBox.Enable = 1;
                app.StandardizeXCheckBox.Value = 0;
                app.StandardizeYCheckBox.Enable = 1;
                app.StandardizeYCheckBox.Value = 0;
                app.FDRThresh.Value = 0.05;
                app.FDRThresh.Enable = 1;
                app.FweThresh.Value = 0.05;
                app.FweThresh.Enable = 1;
                app.UncPThresh.Value = .001;
                app.UncPThresh.Enable = 1;      
                app.ApplyDTLVCCheckBox.Enable = 1;
                app.ApplyDTLVCCheckBox.Value = 0; 
                app.MinFreq.Value = 10;
                app.RegressLesionVolumefromYCheckBox.Enable = 1;
                app.SelectNuisanceRegressorsListBox.Enable = 1;
                app.NuisanceRegressorCheckBox.Enable = 1;                    
                
            end
        end

        % Value changed function: OptimizeHyperparametersCheckBox
        function OptimizeHyperparametersCheckBoxValueChanged(app, event)
            value = app.OptimizeHyperparametersCheckBox.Value;
            if value == 1
                app.HpOptCvType.Enable = "on";
                app.HpOptFolds.Enable = "on";
                app.HpOptRepeats.Enable = "on";
            elseif value == 0
                app.HpOptCvType.Enable = "off";
                app.HpOptFolds.Enable = "off";
                app.HpOptRepeats.Enable = "off";                
            end
        end

        % Value changed function: EnableModelStackingCheckBox
        function EnableModelStackingCheckBoxValueChanged(app, event)
            
        end

        % Button down function: PredictorModalityButtonGroup
        function PredictorModalityButtonGroupButtonDown(app, event)

        end

        % Selection changed function: PredictorModalityButtonGroup
        function PredictorModalityButtonGroupSelectionChanged(app, event)
            if app.MatrixButton.Value == 1
                app.BrainMaskText.Enable = 1;
                app.BrainMaskText.Editable = 1;
                app.BrainMaskText.Placeholder = 'Parcellation table not selected (required to write .node files)';     
                app.SelectBrainMaskButton.Text = 'Select Parcel Table';
                app.SelectBrainMaskButton.Enable = 1;
            elseif app.NIFTIButton.Value == 1
                app.BrainMaskText.Enable = 1;
                app.BrainMaskText.Editable = 1;
                app.BrainMaskText.Placeholder = 'No mask selected';    
                app.SelectBrainMaskButton.Text = 'Select Brain Mask';
                app.SelectBrainMaskButton.Enable = 1;                  
            elseif app.ArrayButton.Value == 1
                app.BrainMaskText.Enable = 0;
                app.BrainMaskText.Editable = 0;
                app.BrainMaskText.Placeholder = 'N/A';
                app.SelectBrainMaskButton.Text = 'N/A';
                app.SelectBrainMaskButton.Enable = 0;     
            elseif app.TableButton.Value == 1
                app.BrainMaskText.Enable = 0;
                app.BrainMaskText.Enable = 0;
                app.BrainMaskText.Editable = 0;
                app.BrainMaskText.Placeholder = 'N/A';
                app.SelectBrainMaskButton.Text = 'N/A';
                app.SelectBrainMaskButton.Enable = 0;                 
            end            
        end

        % Callback function
        function PredictorTypeButtonGroupButtonDown(app, event)
            
        end

        % Value changed function: CVPermCheckBox
        function CVPermCheckBoxValueChanged(app, event)
            value = app.CVPermCheckBox.Value;
            if value == 1
                app.CVPermIter.Enable = 1;
                app.CVPermIter.Value = 50;
            else
                app.CVPermIter.Enable = 0;
                app.CVPermIter.Value = 0;
            end
        end

        % Value changed function: CVType
        function CVTypeValueChanged(app, event)
            value = app.CVType.Value;
            if value == 1
                app.OuterFolds.Enable = 1;
                app.OuterFolds.Value = 5;
                app.HpOptRepeats.Enable = 1;
                app.HpOptRepeats.Value = 20;
            else
                app.OuterFolds.Enable = 0;
                app.OuterFolds.Value = 0;
                app.HpOptRepeats.Enable = 0;
                app.HpOptRepeats.Value = 0;
            end
        end

        % Value changed function: NuisanceRegressorCheckBox
        function NuisanceRegressorCheckBoxValueChanged(app, event)
            value = app.NuisanceRegressorCheckBox.Value;
            if value == 1
                app.SelectNuisanceRegressorsListBox.Enable = 1;
                app.SelectNuisanceRegressorsListBox.Items = setdiff(app.SelectOutcomeVariableYListBox.Items, app.SelectOutcomeVariableYListBox.Value);
            else
                app.SelectNuisanceRegressorsListBox.Enable = 0;
                app.SelectNuisanceRegressorsListBox.Items = string();
                app.SelectNuisanceRegressorsListBox.Value = {};
            end
        end

        % Value changed function: CatYCheckBox
        function CatYCheckBoxValueChanged(app, event)
            value = app.CatYCheckBox.Value;
            if value == 1
                app.RegressLesionVolumefromYCheckBox.Text = 'Regress lesion volume from X';
            else
                app.RegressLesionVolumefromYCheckBox.Text = 'Regress lesion volume from Y';
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 857 886];
            app.UIFigure.Name = 'MATLAB App';

            % Create GeneralModelingOptionsPanel
            app.GeneralModelingOptionsPanel = uipanel(app.UIFigure);
            app.GeneralModelingOptionsPanel.Title = 'General Modeling Options';
            app.GeneralModelingOptionsPanel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.GeneralModelingOptionsPanel.FontWeight = 'bold';
            app.GeneralModelingOptionsPanel.FontSize = 14;
            app.GeneralModelingOptionsPanel.Position = [9 212 833 207];

            % Create StandardizeYCheckBox
            app.StandardizeYCheckBox = uicheckbox(app.GeneralModelingOptionsPanel);
            app.StandardizeYCheckBox.Text = 'Standardize Y';
            app.StandardizeYCheckBox.Position = [6 152 97 22];

            % Create StandardizeXCheckBox
            app.StandardizeXCheckBox = uicheckbox(app.GeneralModelingOptionsPanel);
            app.StandardizeXCheckBox.Text = 'Standardize X';
            app.StandardizeXCheckBox.Position = [108 152 97 22];

            % Create ApplyDTLVCCheckBox
            app.ApplyDTLVCCheckBox = uicheckbox(app.GeneralModelingOptionsPanel);
            app.ApplyDTLVCCheckBox.Text = 'Apply DTLVC';
            app.ApplyDTLVCCheckBox.Position = [219 152 94 22];
            app.ApplyDTLVCCheckBox.Value = true;

            % Create PermApplycFWECheckBox
            app.PermApplycFWECheckBox = uicheckbox(app.GeneralModelingOptionsPanel);
            app.PermApplycFWECheckBox.Text = 'Apply cFWE';
            app.PermApplycFWECheckBox.Position = [413 47 88 22];

            % Create PermVoxPValsCheckBox
            app.PermVoxPValsCheckBox = uicheckbox(app.GeneralModelingOptionsPanel);
            app.PermVoxPValsCheckBox.Enable = 'off';
            app.PermVoxPValsCheckBox.Visible = 'off';
            app.PermVoxPValsCheckBox.Text = 'Voxel p-values';
            app.PermVoxPValsCheckBox.Position = [413 19 100 22];

            % Create BootCIsCheckBox
            app.BootCIsCheckBox = uicheckbox(app.GeneralModelingOptionsPanel);
            app.BootCIsCheckBox.Text = 'Bootstrap CIs';
            app.BootCIsCheckBox.Position = [550 47 95 22];
            app.BootCIsCheckBox.Value = true;

            % Create ConfIntLabel
            app.ConfIntLabel = uilabel(app.GeneralModelingOptionsPanel);
            app.ConfIntLabel.HorizontalAlignment = 'right';
            app.ConfIntLabel.Position = [544 22 73 22];
            app.ConfIntLabel.Text = 'CI Width (%)';

            % Create PermIter
            app.PermIter = uieditfield(app.GeneralModelingOptionsPanel, 'numeric');
            app.PermIter.Limits = [0 Inf];
            app.PermIter.Position = [492 71 46 22];
            app.PermIter.Value = 1000;

            % Create ConfIntWidth
            app.ConfIntWidth = uieditfield(app.GeneralModelingOptionsPanel, 'numeric');
            app.ConfIntWidth.Limits = [0 Inf];
            app.ConfIntWidth.Position = [625 22 37 22];
            app.ConfIntWidth.Value = 95;

            % Create StratifyTrainTestSplitsonYCheckBox
            app.StratifyTrainTestSplitsonYCheckBox = uicheckbox(app.GeneralModelingOptionsPanel);
            app.StratifyTrainTestSplitsonYCheckBox.ValueChangedFcn = createCallbackFcn(app, @StratifyTrainTestSplitsonYCheckBoxValueChanged, true);
            app.StratifyTrainTestSplitsonYCheckBox.Text = 'Stratify Train/Test Splits on Y';
            app.StratifyTrainTestSplitsonYCheckBox.Position = [322 152 175 22];

            % Create PermCheckBox
            app.PermCheckBox = uicheckbox(app.GeneralModelingOptionsPanel);
            app.PermCheckBox.ValueChangedFcn = createCallbackFcn(app, @PermCheckBoxValueChanged, true);
            app.PermCheckBox.Text = 'Run Permutations';
            app.PermCheckBox.Position = [414 94 118 22];
            app.PermCheckBox.Value = true;

            % Create BootCheckBox
            app.BootCheckBox = uicheckbox(app.GeneralModelingOptionsPanel);
            app.BootCheckBox.ValueChangedFcn = createCallbackFcn(app, @BootCheckBoxValueChanged, true);
            app.BootCheckBox.Text = 'Run Bootstraps';
            app.BootCheckBox.Position = [550 94 105 22];
            app.BootCheckBox.Value = true;

            % Create PermutationsEditFieldLabel
            app.PermutationsEditFieldLabel = uilabel(app.GeneralModelingOptionsPanel);
            app.PermutationsEditFieldLabel.HorizontalAlignment = 'right';
            app.PermutationsEditFieldLabel.Position = [408 71 76 22];
            app.PermutationsEditFieldLabel.Text = 'Permutations';

            % Create BootstrapLabel
            app.BootstrapLabel = uilabel(app.GeneralModelingOptionsPanel);
            app.BootstrapLabel.HorizontalAlignment = 'right';
            app.BootstrapLabel.Position = [545 71 62 22];
            app.BootstrapLabel.Text = 'Bootstraps';

            % Create BootAnalysesPanel
            app.BootAnalysesPanel = uipanel(app.GeneralModelingOptionsPanel);
            app.BootAnalysesPanel.BorderType = 'none';
            app.BootAnalysesPanel.Title = 'Bootstrap Analyses';
            app.BootAnalysesPanel.FontWeight = 'bold';
            app.BootAnalysesPanel.Position = [548 113 135 30];

            % Create BootIter
            app.BootIter = uieditfield(app.GeneralModelingOptionsPanel, 'numeric');
            app.BootIter.Limits = [0 Inf];
            app.BootIter.Position = [616 71 46 22];
            app.BootIter.Value = 1000;

            % Create CrossValidationOptionsNestedifHyperParameterTuningEnabledPanel
            app.CrossValidationOptionsNestedifHyperParameterTuningEnabledPanel = uipanel(app.UIFigure);
            app.CrossValidationOptionsNestedifHyperParameterTuningEnabledPanel.Title = 'Cross-Validation Options (Nested if Hyper-Parameter Tuning Enabled)';
            app.CrossValidationOptionsNestedifHyperParameterTuningEnabledPanel.FontWeight = 'bold';
            app.CrossValidationOptionsNestedifHyperParameterTuningEnabledPanel.FontSize = 14;
            app.CrossValidationOptionsNestedifHyperParameterTuningEnabledPanel.Position = [10 9 544 174];

            % Create CrossValidationSettingsPanel
            app.CrossValidationSettingsPanel = uipanel(app.CrossValidationOptionsNestedifHyperParameterTuningEnabledPanel);
            app.CrossValidationSettingsPanel.BorderType = 'none';
            app.CrossValidationSettingsPanel.Title = 'Cross-Validation Settings';
            app.CrossValidationSettingsPanel.FontWeight = 'bold';
            app.CrossValidationSettingsPanel.Position = [182 28 205 120];

            % Create CVStrategyLabel_2
            app.CVStrategyLabel_2 = uilabel(app.CrossValidationSettingsPanel);
            app.CVStrategyLabel_2.Position = [11 70 70 22];
            app.CVStrategyLabel_2.Text = 'CV Strategy';

            % Create NumberofOuterFoldsLabel
            app.NumberofOuterFoldsLabel = uilabel(app.CrossValidationSettingsPanel);
            app.NumberofOuterFoldsLabel.Position = [11 41 127 22];
            app.NumberofOuterFoldsLabel.Text = 'Number of Outer Folds';

            % Create NumberofRepeatsLabel_2
            app.NumberofRepeatsLabel_2 = uilabel(app.CrossValidationSettingsPanel);
            app.NumberofRepeatsLabel_2.Position = [11 11 109 22];
            app.NumberofRepeatsLabel_2.Text = 'Number of Repeats';

            % Create CVType
            app.CVType = uidropdown(app.CrossValidationSettingsPanel);
            app.CVType.Items = {'KFold', 'LOOCV'};
            app.CVType.Editable = 'on';
            app.CVType.ValueChangedFcn = createCallbackFcn(app, @CVTypeValueChanged, true);
            app.CVType.BackgroundColor = [1 1 1];
            app.CVType.Position = [131 70 69 22];
            app.CVType.Value = 'KFold';

            % Create OuterFolds
            app.OuterFolds = uieditfield(app.CrossValidationSettingsPanel, 'numeric');
            app.OuterFolds.Position = [167 41 33 22];
            app.OuterFolds.Value = 5;

            % Create OuterRepeats
            app.OuterRepeats = uieditfield(app.CrossValidationSettingsPanel, 'numeric');
            app.OuterRepeats.Position = [167 10 33 22];
            app.OuterRepeats.Value = 5;

            % Create RunCrossValidationCheckBox
            app.RunCrossValidationCheckBox = uicheckbox(app.CrossValidationOptionsNestedifHyperParameterTuningEnabledPanel);
            app.RunCrossValidationCheckBox.ValueChangedFcn = createCallbackFcn(app, @RunCrossValidationCheckBoxValueChanged, true);
            app.RunCrossValidationCheckBox.Text = 'Run Cross-Validation';
            app.RunCrossValidationCheckBox.Position = [7 124 135 22];
            app.RunCrossValidationCheckBox.Value = true;

            % Create CVPermPanel
            app.CVPermPanel = uipanel(app.CrossValidationOptionsNestedifHyperParameterTuningEnabledPanel);
            app.CVPermPanel.BorderType = 'none';
            app.CVPermPanel.Title = 'Permutation Testing';
            app.CVPermPanel.FontWeight = 'bold';
            app.CVPermPanel.Position = [395 89 143 59];

            % Create CVPermCheckBox
            app.CVPermCheckBox = uicheckbox(app.CVPermPanel);
            app.CVPermCheckBox.ValueChangedFcn = createCallbackFcn(app, @CVPermCheckBoxValueChanged, true);
            app.CVPermCheckBox.Text = 'Run Permutations';
            app.CVPermCheckBox.Position = [9 9 118 22];
            app.CVPermCheckBox.Value = true;

            % Create PermutationsEditField_2Label
            app.PermutationsEditField_2Label = uilabel(app.CrossValidationOptionsNestedifHyperParameterTuningEnabledPanel);
            app.PermutationsEditField_2Label.HorizontalAlignment = 'right';
            app.PermutationsEditField_2Label.Position = [397 68 76 22];
            app.PermutationsEditField_2Label.Text = 'Permutations';

            % Create CVPermIter
            app.CVPermIter = uieditfield(app.CrossValidationOptionsNestedifHyperParameterTuningEnabledPanel, 'numeric');
            app.CVPermIter.Position = [490 68 37 22];
            app.CVPermIter.Value = 50;

            % Create EnableModelStackingCheckBox
            app.EnableModelStackingCheckBox = uicheckbox(app.CrossValidationOptionsNestedifHyperParameterTuningEnabledPanel);
            app.EnableModelStackingCheckBox.ValueChangedFcn = createCallbackFcn(app, @EnableModelStackingCheckBoxValueChanged, true);
            app.EnableModelStackingCheckBox.Text = 'Enable Model Stacking';
            app.EnableModelStackingCheckBox.Position = [7 97 145 22];
            app.EnableModelStackingCheckBox.Value = true;

            % Create PredictorModalityButtonGroup
            app.PredictorModalityButtonGroup = uibuttongroup(app.UIFigure);
            app.PredictorModalityButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @PredictorModalityButtonGroupSelectionChanged, true);
            app.PredictorModalityButtonGroup.TitlePosition = 'centertop';
            app.PredictorModalityButtonGroup.Title = 'Predictor Modality';
            app.PredictorModalityButtonGroup.ButtonDownFcn = createCallbackFcn(app, @PredictorModalityButtonGroupButtonDown, true);
            app.PredictorModalityButtonGroup.FontWeight = 'bold';
            app.PredictorModalityButtonGroup.Position = [10 766 125 69];

            % Create NIFTIButton
            app.NIFTIButton = uiradiobutton(app.PredictorModalityButtonGroup);
            app.NIFTIButton.Text = 'NIFTI';
            app.NIFTIButton.Position = [11 23 52 22];
            app.NIFTIButton.Value = true;

            % Create MatrixButton
            app.MatrixButton = uiradiobutton(app.PredictorModalityButtonGroup);
            app.MatrixButton.Text = 'Matrix';
            app.MatrixButton.Position = [65 22 65 22];

            % Create ArrayButton
            app.ArrayButton = uiradiobutton(app.PredictorModalityButtonGroup);
            app.ArrayButton.Text = 'Array';
            app.ArrayButton.Position = [11 1 65 22];

            % Create SelectBrainMaskButton
            app.SelectBrainMaskButton = uibutton(app.UIFigure, 'push');
            app.SelectBrainMaskButton.ButtonPushedFcn = createCallbackFcn(app, @SelectBrainMaskButtonPushed, true);
            app.SelectBrainMaskButton.BackgroundColor = [0.651 0.651 0.651];
            app.SelectBrainMaskButton.Position = [163 700 135 23];
            app.SelectBrainMaskButton.Text = 'Select Brain Mask';

            % Create HyperParameterTuningPanel
            app.HyperParameterTuningPanel = uipanel(app.UIFigure);
            app.HyperParameterTuningPanel.BorderType = 'none';
            app.HyperParameterTuningPanel.Title = 'Hyper-Parameter Tuning';
            app.HyperParameterTuningPanel.FontWeight = 'bold';
            app.HyperParameterTuningPanel.Position = [199 221 205 134];

            % Create HpOptCVStrategyLabel
            app.HpOptCVStrategyLabel = uilabel(app.HyperParameterTuningPanel);
            app.HpOptCVStrategyLabel.Position = [11 65 70 22];
            app.HpOptCVStrategyLabel.Text = 'CV Strategy';

            % Create HpOptNumberofFoldsLabel
            app.HpOptNumberofFoldsLabel = uilabel(app.HyperParameterTuningPanel);
            app.HpOptNumberofFoldsLabel.Position = [11 36 94 22];
            app.HpOptNumberofFoldsLabel.Text = 'Number of Folds';

            % Create HpOptNumberofRepeatsLabel
            app.HpOptNumberofRepeatsLabel = uilabel(app.HyperParameterTuningPanel);
            app.HpOptNumberofRepeatsLabel.Position = [11 9 109 22];
            app.HpOptNumberofRepeatsLabel.Text = 'Number of Repeats';

            % Create HpOptCvType
            app.HpOptCvType = uidropdown(app.HyperParameterTuningPanel);
            app.HpOptCvType.Items = {'KFold', 'LOOCV'};
            app.HpOptCvType.ValueChangedFcn = createCallbackFcn(app, @HpOptCvTypeValueChanged, true);
            app.HpOptCvType.Position = [131 65 69 22];
            app.HpOptCvType.Value = 'KFold';

            % Create HpOptFolds
            app.HpOptFolds = uieditfield(app.HyperParameterTuningPanel, 'numeric');
            app.HpOptFolds.Position = [160 36 40 22];
            app.HpOptFolds.Value = 5;

            % Create HpOptRepeats
            app.HpOptRepeats = uieditfield(app.HyperParameterTuningPanel, 'numeric');
            app.HpOptRepeats.Position = [160 8 40 22];
            app.HpOptRepeats.Value = 5;

            % Create OptimizeHyperparametersCheckBox
            app.OptimizeHyperparametersCheckBox = uicheckbox(app.HyperParameterTuningPanel);
            app.OptimizeHyperparametersCheckBox.ValueChangedFcn = createCallbackFcn(app, @OptimizeHyperparametersCheckBoxValueChanged, true);
            app.OptimizeHyperparametersCheckBox.Text = 'Optimize Hyper-parameters';
            app.OptimizeHyperparametersCheckBox.Position = [11 89 169 22];
            app.OptimizeHyperparametersCheckBox.Value = true;

            % Create SelectDataDirectoryButton
            app.SelectDataDirectoryButton = uibutton(app.UIFigure, 'push');
            app.SelectDataDirectoryButton.ButtonPushedFcn = createCallbackFcn(app, @SelectDataDirectoryButtonPushed, true);
            app.SelectDataDirectoryButton.BackgroundColor = [0.651 0.651 0.651];
            app.SelectDataDirectoryButton.Position = [163 743 136 24];
            app.SelectDataDirectoryButton.Text = 'Select Data Directory';

            % Create SelectCSVfileButton
            app.SelectCSVfileButton = uibutton(app.UIFigure, 'push');
            app.SelectCSVfileButton.ButtonPushedFcn = createCallbackFcn(app, @SelectCSVfileButtonPushed, true);
            app.SelectCSVfileButton.BackgroundColor = [0.651 0.651 0.651];
            app.SelectCSVfileButton.Position = [164 811 134 23];
            app.SelectCSVfileButton.Text = 'Select CSV file';

            % Create MinimumDataThresholdsPanel
            app.MinimumDataThresholdsPanel = uipanel(app.UIFigure);
            app.MinimumDataThresholdsPanel.BorderType = 'none';
            app.MinimumDataThresholdsPanel.Title = 'Minimum Data Thresholds';
            app.MinimumDataThresholdsPanel.FontWeight = 'bold';
            app.MinimumDataThresholdsPanel.Position = [12 266 180 89];

            % Create MinimumFrequencyLabel
            app.MinimumFrequencyLabel = uilabel(app.MinimumDataThresholdsPanel);
            app.MinimumFrequencyLabel.Position = [10 10 114 22];
            app.MinimumFrequencyLabel.Text = 'Minimum Frequency';

            % Create MinFreq
            app.MinFreq = uieditfield(app.MinimumDataThresholdsPanel, 'numeric');
            app.MinFreq.Position = [142 10 39 22];
            app.MinFreq.Value = 10;

            % Create MinimumObservationLabel
            app.MinimumObservationLabel = uilabel(app.MinimumDataThresholdsPanel);
            app.MinimumObservationLabel.Position = [10 41 122 22];
            app.MinimumObservationLabel.Text = 'Minimum Observation';

            % Create MinObs
            app.MinObs = uieditfield(app.MinimumDataThresholdsPanel, 'numeric');
            app.MinObs.Position = [142 41 39 22];
            app.MinObs.Value = 1;

            % Create PredictorTypeButtonGroup
            app.PredictorTypeButtonGroup = uibuttongroup(app.UIFigure);
            app.PredictorTypeButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @PredictorTypeButtonGroupSelectionChanged, true);
            app.PredictorTypeButtonGroup.TitlePosition = 'centertop';
            app.PredictorTypeButtonGroup.Title = 'Predictor Type';
            app.PredictorTypeButtonGroup.FontWeight = 'bold';
            app.PredictorTypeButtonGroup.Position = [10 701 124 51];

            % Create LesionButton
            app.LesionButton = uiradiobutton(app.PredictorTypeButtonGroup);
            app.LesionButton.Text = 'Lesion';
            app.LesionButton.Position = [11 5 58 22];
            app.LesionButton.Value = true;

            % Create OtherPredButton
            app.OtherPredButton = uiradiobutton(app.PredictorTypeButtonGroup);
            app.OtherPredButton.Text = 'Other';
            app.OtherPredButton.Position = [69 5 65 22];

            % Create CSVTtext
            app.CSVTtext = uitextarea(app.UIFigure);
            app.CSVTtext.ValueChangedFcn = createCallbackFcn(app, @SelectCSVfileButtonPushed, true);
            app.CSVTtext.Placeholder = 'No file selected';
            app.CSVTtext.Position = [308 811 533 24];

            % Create ImDirText
            app.ImDirText = uitextarea(app.UIFigure);
            app.ImDirText.Placeholder = 'No directory selected';
            app.ImDirText.Position = [308 743 533 24];

            % Create BrainMaskText
            app.BrainMaskText = uitextarea(app.UIFigure);
            app.BrainMaskText.Placeholder = 'No mask selected (NIFTI Modality Only)';
            app.BrainMaskText.Position = [308 699 533 25];

            % Create SelectModelLabel
            app.SelectModelLabel = uilabel(app.UIFigure);
            app.SelectModelLabel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.SelectModelLabel.HorizontalAlignment = 'right';
            app.SelectModelLabel.FontSize = 14;
            app.SelectModelLabel.Position = [9 446 90 22];
            app.SelectModelLabel.Text = 'Select Model:';

            % Create SelectModelDropDown
            app.SelectModelDropDown = uidropdown(app.UIFigure);
            app.SelectModelDropDown.Items = {'plsr', 'pls_da', 'ridge', 'logistic_ridge', 'linsvr', 'linsvc', 'kernsvr', 'kernsvc', 'rensemble', 'censemble', 'municorr', 'ttest', 'bmunz', 'munilr', 'muniolsr', 'prop_sub'};
            app.SelectModelDropDown.ValueChangedFcn = createCallbackFcn(app, @SelectModelDropDownValueChanged, true);
            app.SelectModelDropDown.FontSize = 14;
            app.SelectModelDropDown.BackgroundColor = [0.651 0.651 0.651];
            app.SelectModelDropDown.Position = [103 442 100 30];
            app.SelectModelDropDown.Value = 'plsr';

            % Create ResultOptionsPanel
            app.ResultOptionsPanel = uipanel(app.UIFigure);
            app.ResultOptionsPanel.Title = 'Result Options';
            app.ResultOptionsPanel.FontWeight = 'bold';
            app.ResultOptionsPanel.FontSize = 14;
            app.ResultOptionsPanel.Position = [595 80 204 77];

            % Create WriteBootstrapImagesCheckBox
            app.WriteBootstrapImagesCheckBox = uicheckbox(app.ResultOptionsPanel);
            app.WriteBootstrapImagesCheckBox.Text = 'Write Bootstrap Images';
            app.WriteBootstrapImagesCheckBox.Position = [6 3 147 22];
            app.WriteBootstrapImagesCheckBox.Value = true;

            % Create WritePermutationImagesCheckBox
            app.WritePermutationImagesCheckBox = uicheckbox(app.ResultOptionsPanel);
            app.WritePermutationImagesCheckBox.Text = 'Write Permutation Images';
            app.WritePermutationImagesCheckBox.Position = [6 27 161 22];
            app.WritePermutationImagesCheckBox.Value = true;

            % Create FitExplanatoryModelCheckBox
            app.FitExplanatoryModelCheckBox = uicheckbox(app.UIFigure);
            app.FitExplanatoryModelCheckBox.Text = 'Fit Explanatory Model';
            app.FitExplanatoryModelCheckBox.FontSize = 14;
            app.FitExplanatoryModelCheckBox.Position = [218 444 157 22];
            app.FitExplanatoryModelCheckBox.Value = true;

            % Create UseParallelProcessingRecommendedCheckBox
            app.UseParallelProcessingRecommendedCheckBox = uicheckbox(app.UIFigure);
            app.UseParallelProcessingRecommendedCheckBox.Text = 'Use Parallel Processing (Recommended)';
            app.UseParallelProcessingRecommendedCheckBox.FontSize = 14;
            app.UseParallelProcessingRecommendedCheckBox.Position = [386 444 282 22];
            app.UseParallelProcessingRecommendedCheckBox.Value = true;

            % Create MultivariateModelingToolAnalysisConfigurationLabel
            app.MultivariateModelingToolAnalysisConfigurationLabel = uilabel(app.UIFigure);
            app.MultivariateModelingToolAnalysisConfigurationLabel.FontSize = 16;
            app.MultivariateModelingToolAnalysisConfigurationLabel.FontWeight = 'bold';
            app.MultivariateModelingToolAnalysisConfigurationLabel.Position = [11 860 503 22];
            app.MultivariateModelingToolAnalysisConfigurationLabel.Text = 'Multivariate Modeling Tool Analysis Configuration ';

            % Create RunAnalysisButton
            app.RunAnalysisButton = uibutton(app.UIFigure, 'push');
            app.RunAnalysisButton.ButtonPushedFcn = createCallbackFcn(app, @RunAnalysisButtonPushed, true);
            app.RunAnalysisButton.BackgroundColor = [0.651 0.651 0.651];
            app.RunAnalysisButton.FontSize = 18;
            app.RunAnalysisButton.FontWeight = 'bold';
            app.RunAnalysisButton.Position = [588 10 224 48];
            app.RunAnalysisButton.Text = 'Run Analysis';

            % Create SelectOutputDirectoryButton
            app.SelectOutputDirectoryButton = uibutton(app.UIFigure, 'push');
            app.SelectOutputDirectoryButton.ButtonPushedFcn = createCallbackFcn(app, @SelectOutputDirectoryButtonPushed, true);
            app.SelectOutputDirectoryButton.BackgroundColor = [0.651 0.651 0.651];
            app.SelectOutputDirectoryButton.Position = [163 656 135 24];
            app.SelectOutputDirectoryButton.Text = 'Select Output Directory';

            % Create OutDirText
            app.OutDirText = uitextarea(app.UIFigure);
            app.OutDirText.Placeholder = 'No directory selected';
            app.OutDirText.Position = [308 657 533 23];

            % Create RegistryFlag
            app.RegistryFlag = uicheckbox(app.UIFigure);
            app.RegistryFlag.Text = 'Registry Patients (select if subject IDs are provided as 4-digit numeric codes in CSV file)';
            app.RegistryFlag.FontSize = 14;
            app.RegistryFlag.Position = [167 782 580 22];

            % Create PermutationTestingPanel
            app.PermutationTestingPanel = uipanel(app.UIFigure);
            app.PermutationTestingPanel.BorderType = 'none';
            app.PermutationTestingPanel.Title = 'Permutation Testing';
            app.PermutationTestingPanel.FontWeight = 'bold';
            app.PermutationTestingPanel.Position = [413 325 135 30];

            % Create StatisticalThresholdsPanel
            app.StatisticalThresholdsPanel = uipanel(app.UIFigure);
            app.StatisticalThresholdsPanel.BorderType = 'none';
            app.StatisticalThresholdsPanel.Title = 'Statistical Thresholds';
            app.StatisticalThresholdsPanel.FontWeight = 'bold';
            app.StatisticalThresholdsPanel.Position = [699 231 135 124];

            % Create FWEEditFieldLabel
            app.FWEEditFieldLabel = uilabel(app.StatisticalThresholdsPanel);
            app.FWEEditFieldLabel.HorizontalAlignment = 'right';
            app.FWEEditFieldLabel.Position = [1 39 32 22];
            app.FWEEditFieldLabel.Text = 'FWE';

            % Create FweThresh
            app.FweThresh = uieditfield(app.StatisticalThresholdsPanel, 'numeric');
            app.FweThresh.Limits = [0 1];
            app.FweThresh.Position = [39 39 44 22];
            app.FweThresh.Value = 0.05;

            % Create UncEditFieldLabel
            app.UncEditFieldLabel = uilabel(app.StatisticalThresholdsPanel);
            app.UncEditFieldLabel.HorizontalAlignment = 'right';
            app.UncEditFieldLabel.Position = [3 74 30 22];
            app.UncEditFieldLabel.Text = 'Unc.';

            % Create UncPThresh
            app.UncPThresh = uieditfield(app.StatisticalThresholdsPanel, 'numeric');
            app.UncPThresh.Limits = [0 1];
            app.UncPThresh.Position = [39 74 44 22];
            app.UncPThresh.Value = 0.001;

            % Create FDRThresh
            app.FDRThresh = uieditfield(app.StatisticalThresholdsPanel, 'numeric');
            app.FDRThresh.Limits = [0 1];
            app.FDRThresh.Position = [39 6 44 22];
            app.FDRThresh.Value = 0.05;

            % Create FDRLabel
            app.FDRLabel = uilabel(app.StatisticalThresholdsPanel);
            app.FDRLabel.HorizontalAlignment = 'right';
            app.FDRLabel.Position = [1 8 30 22];
            app.FDRLabel.Text = 'FDR';

            % Create MisclassificationCostsPanel
            app.MisclassificationCostsPanel = uipanel(app.UIFigure);
            app.MisclassificationCostsPanel.BorderType = 'none';
            app.MisclassificationCostsPanel.Title = 'Misclassification Costs';
            app.MisclassificationCostsPanel.FontWeight = 'bold';
            app.MisclassificationCostsPanel.Position = [15 10 168 90];

            % Create Group1EditField_4Label
            app.Group1EditField_4Label = uilabel(app.MisclassificationCostsPanel);
            app.Group1EditField_4Label.HorizontalAlignment = 'right';
            app.Group1EditField_4Label.Position = [15 13 52 22];
            app.Group1EditField_4Label.Text = 'Group -1';

            % Create MisClasCostG2
            app.MisClasCostG2 = uieditfield(app.MisclassificationCostsPanel, 'numeric');
            app.MisClasCostG2.Limits = [0 100];
            app.MisClasCostG2.Enable = 'off';
            app.MisClasCostG2.Position = [74 13 22 22];
            app.MisClasCostG2.Value = 1;

            % Create Group1EditField_3Label
            app.Group1EditField_3Label = uilabel(app.MisclassificationCostsPanel);
            app.Group1EditField_3Label.HorizontalAlignment = 'right';
            app.Group1EditField_3Label.Position = [13 40 55 22];
            app.Group1EditField_3Label.Text = 'Group +1';

            % Create MisClasCostG1
            app.MisClasCostG1 = uieditfield(app.MisclassificationCostsPanel, 'numeric');
            app.MisClasCostG1.Limits = [0 1000];
            app.MisClasCostG1.Enable = 'off';
            app.MisClasCostG1.Position = [74 40 22 22];
            app.MisClasCostG1.Value = 1;

            % Create SelectOutcomeVariableYListBoxLabel
            app.SelectOutcomeVariableYListBoxLabel = uilabel(app.UIFigure);
            app.SelectOutcomeVariableYListBoxLabel.HorizontalAlignment = 'right';
            app.SelectOutcomeVariableYListBoxLabel.Position = [50 603 156 22];
            app.SelectOutcomeVariableYListBoxLabel.Text = 'Select Outcome Variable (Y)';

            % Create SelectOutcomeVariableYListBox
            app.SelectOutcomeVariableYListBox = uilistbox(app.UIFigure);
            app.SelectOutcomeVariableYListBox.Items = {};
            app.SelectOutcomeVariableYListBox.Position = [221 585 620 56];
            app.SelectOutcomeVariableYListBox.Value = {};

            % Create SelectNuisanceRegressorsListBoxLabel
            app.SelectNuisanceRegressorsListBoxLabel = uilabel(app.UIFigure);
            app.SelectNuisanceRegressorsListBoxLabel.HorizontalAlignment = 'right';
            app.SelectNuisanceRegressorsListBoxLabel.Position = [34 533 156 22];
            app.SelectNuisanceRegressorsListBoxLabel.Text = 'Select Nuisance Regressors';

            % Create SelectNuisanceRegressorsListBox
            app.SelectNuisanceRegressorsListBox = uilistbox(app.UIFigure);
            app.SelectNuisanceRegressorsListBox.Items = {};
            app.SelectNuisanceRegressorsListBox.Multiselect = 'on';
            app.SelectNuisanceRegressorsListBox.Position = [221 514 620 56];
            app.SelectNuisanceRegressorsListBox.Value = {};

            % Create RegressLesionVolumefromYCheckBox
            app.RegressLesionVolumefromYCheckBox = uicheckbox(app.UIFigure);
            app.RegressLesionVolumefromYCheckBox.Text = 'Regress Lesion Volume from Y';
            app.RegressLesionVolumefromYCheckBox.Position = [221 489 187 22];

            % Create NuisanceRegressorCheckBox
            app.NuisanceRegressorCheckBox = uicheckbox(app.UIFigure);
            app.NuisanceRegressorCheckBox.ValueChangedFcn = createCallbackFcn(app, @NuisanceRegressorCheckBoxValueChanged, true);
            app.NuisanceRegressorCheckBox.Text = '';
            app.NuisanceRegressorCheckBox.Position = [196 532 25 22];

            % Create CatYCheckBox
            app.CatYCheckBox = uicheckbox(app.UIFigure);
            app.CatYCheckBox.ValueChangedFcn = createCallbackFcn(app, @CatYCheckBoxValueChanged, true);
            app.CatYCheckBox.Enable = 'off';
            app.CatYCheckBox.Visible = 'off';
            app.CatYCheckBox.Text = 'Y is Grouping Variable';
            app.CatYCheckBox.Position = [217 421 141 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = run_modeling_gui

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end