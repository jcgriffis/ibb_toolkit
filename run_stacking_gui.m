classdef run_stacking_gui < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        PlottingOptionsPanel          matlab.ui.container.Panel
        ComparisonPlotsPanel          matlab.ui.container.Panel
        MakePlots                     matlab.ui.control.CheckBox
        GeneralModelingOptionsPanel   matlab.ui.container.Panel
        HyperParameterTuningPanel     matlab.ui.container.Panel
        OptimizeHyperparametersCheckBox  matlab.ui.control.CheckBox
        Repeats                       matlab.ui.control.NumericEditField
        Folds                         matlab.ui.control.NumericEditField
        CVType                        matlab.ui.control.DropDown
        NRepeatsLabel                 matlab.ui.control.Label
        NFoldsLabel                   matlab.ui.control.Label
        CVStrategyLabel               matlab.ui.control.Label
        StandardizeXCheckBox          matlab.ui.control.CheckBox
        StandardizeYCheckBox          matlab.ui.control.CheckBox
        CrossValidationOptionsPanel   matlab.ui.container.Panel
        MisclassificationCostsLabel   matlab.ui.control.Label
        NegGroupCost                  matlab.ui.control.NumericEditField
        Group1EditFieldLabel          matlab.ui.control.Label
        PosGroupCost                  matlab.ui.control.NumericEditField
        Group1EditField_2Label        matlab.ui.control.Label
        PermIter                      matlab.ui.control.NumericEditField
        PermutationsEditField_2Label  matlab.ui.control.Label
        PermutationTestingPanel       matlab.ui.container.Panel
        RunPerm                       matlab.ui.control.CheckBox
        FileSelText                   matlab.ui.control.TextArea
        OutDirText                    matlab.ui.control.TextArea
        SelectOutputDirectoryButton   matlab.ui.control.Button
        LoadTextFileButton            matlab.ui.control.Button
        RunAnalysisButton             matlab.ui.control.Button
        StackedModelingAnalysisConfigurationLabel  matlab.ui.control.Label
        UseParallelProcessingRecommendedCheckBox  matlab.ui.control.CheckBox
        SelectModelDropDown           matlab.ui.control.DropDown
        SelectModelLabel              matlab.ui.control.Label
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Callback function
        function SelectCSVfileButtonPushed(app, event)

        end

        % Callback function
        function SelectDataDirectoryButtonPushed(app, event)

        end

        % Callback function
        function SelectBrainMaskButtonPushed(app, event)


        end

        % Button pushed function: RunAnalysisButton
        function RunAnalysisButtonPushed(app, event)
            
            % Get model directory
            model_dir = fileparts(which('run_stacking_gui.mlapp'));
            cfg.model_dir = fullfile(model_dir, 'model_code');
            
            % Add core functions to path
            addpath(fullfile(cfg.model_dir, 'common'));
            addpath(fullfile(cfg.model_dir, 'stacked_models'));           

            % Set up basic cfg file
            cfg.model_spec = string(app.SelectModelDropDown.Value);
            cfg.optimize_hyperparams = app.OptimizeHyperparametersCheckBox.Value;
            cfg.cross_validation = 1;      
            cfg.out_dir = string(app.OutDirText.Value);
            cfg.modality = 'other';

            % Get default cfg
            cfg = get_default_stack_cfg(cfg);

            %%%%%% Update defaults with user choices

            %%%% General analysis params

            % Standardization
            if app.StandardizeXCheckBox.Value == 1
                cfg.standardize = 1;
            elseif app.StandardizeYCheckBox.Value == 1 && app.StandardizeXCheckBox.Value == 0
                cfg.standardize = 2;
            elseif app.StandardizeYCheckBox.Value == 1 && app.StandardizeXCheckBox.Value == 1
                cfg.standardize = 3;
            else
                cfg.standardize = 0;
            end            

            % Hyperparameter optimization
            cfg.hp_opt.cv_type = string(app.CVType.Value);        
            cfg.hp_opt.folds = app.Folds.Value;
            cfg.hp_opt.repeats = app.Repeats.Value;           
            
            % Cross-validation
            cfg.cv.stacked_model = 1;
            cfg.cv.permutation = app.RunPerm.Value;
            cfg.cv.n_perm = app.PermIter.Value;
            cfg.cv.save_stacked_results = 1;
            
            % Classification Costs
            if cfg.cat_Y == 1
                cfg.cost(1,2) = app.PosGroupCost.Value;
                cfg.cost(2,1) = app.NegGroupCost.Value;
            end

            % Parallel processing
            cfg.parallel = app.UseParallelProcessingRecommendedCheckBox.Value;
            
            % Load and format base model data
            cfg.model_paths = readlines(app.FileSelText.Value{1})';
       
            % Run stacked modeling
            [stacked_results, base_results] = run_model_stacking(cfg);
            
            % Save results
            cd(cfg.out_dir);
            save('stacked_results.mat', 'stacked_results', 'base_results');

            % Make plots if indicated
            if app.MakePlots.Value == 1
                if cfg.cat_Y == 0
                    % R-squared plot
                    subplot(1,2,1);
                    r2_data = [stacked_results.avg.r2_ss, cell2mat(base_results.avg.r2_ss)];
                    boxplot(r2_data);
                    xtl = cell(1,size(r2_data,2));
                    for i = 1:length(r2_data)
                        if i == 1
                            xtl{i} = 'Stacked';
                        else
                            xtl{i} = ['Model ' num2str(i-1)];
                        end
                    end
                    set(gca, 'XtickLabel', xtl);
                    ylabel('Mean Out-of-Fold R-Squared');
                    title('CV R^2');
                    grid on;
                    % Correlation plot
                    subplot(1,2,2);
                    boxplot([stacked_results.avg.corr cell2mat(base_results.avg.corr)]);
                    set(gca, 'XtickLabel', xtl);
                    ylabel('Mean Out-of-Fold Correlation');
                    title('CV Correlation');
                    grid on;
                else
                    % AUC plot
                    subplot(1,3,1);
                    auc_data = [stacked_results.avg.roc_auc, cell2mat(base_results.avg.roc_auc)];
                    boxplot(auc_data);
                    xtl = cell(1,size(auc_data,2));
                    base_cr1 = zeros(size(base_results.avg.classrate{1}(:,1),1), length(base_results.avg.classrate));
                    base_cr2 = zeros(size(base_results.avg.classrate{1}(:,1),1), length(base_results.avg.classrate));                    
                    for i = 1:size(auc_data,2)
                        if i == 1
                            xtl{i} = 'Stacked';
                        else
                            xtl{i} = ['Model ' num2str(i-1)];
                            base_cr1(:,i-1) = base_results.avg.classrate{i-1}(:,1);
                            base_cr2(:,i-1) = base_results.avg.classrate{i-1}(:,2);                            
                        end
                    end
                    set(gca, 'XtickLabel', xtl);
                    ylabel('Mean Out-of-Fold AUC');
                    title('CV AUC');
                    grid on;
                    % Classification plot
                    subplot(1,3,2);
                    boxplot([stacked_results.avg.classrate(:,1) base_cr1]);
                    set(gca, 'XtickLabel', xtl);
                    ylabel('Mean Out-of-Fold Classification Rate');
                    title('CV Accuracy (Group -1)');
                    grid on;
                    % Classification plot
                    subplot(1,3,3);
                    boxplot([stacked_results.avg.classrate(:,2) base_cr2]);
                    set(gca, 'XtickLabel', xtl);
                    ylabel('Mean Out-of-Fold Classification Rate');
                    title('CV Accuracy (Group +1)');
                    grid on;            
                end
            end
        end

        % Callback function
        function SelectOutputDirectoryButtonPushed(app, event)
             [out_dir] = uigetdir(pwd, 'Select directory to save results');
            if out_dir == 0
                error('No directory selected; please try again');
            else
                app.OutDirText.Value = out_dir;
            end           
        end

        % Callback function
        function cv_typeValueChanged(app, event)

        end

        % Callback function
        function StratifyTrainTestSplitsonYCheckBoxValueChanged(app, event)

        end

        % Callback function
        function RunPermutationsCheckBoxValueChanged(app, event)

        end

        % Callback function
        function RunBootstrapsCheckBoxValueChanged(app, event)
           
        end

        % Callback function
        function PredictorTypeButtonGroupSelectionChanged(app, event)

        end

        % Callback function
        function RunCrossValidationCheckBoxValueChanged(app, event)

        end

        % Value changed function: SelectModelDropDown
        function SelectModelDropDownValueChanged(app, event)

            if ~contains(app.SelectModelDropDown.Value, {'mean', 'mean_score'})
                app.Repeats.Value = 5;
                app.Repeats.Enable = 5;
            end
       
        end

        % Callback function
        function OptimizeHyperparametersCheckBoxValueChanged(app, event)
            value = app.OptimizeHyperparametersCheckBox.Value;
            if value == 1
                app.CVType.Enable = "on";
                app.Folds.Enable = "on";
                app.Repeats.Enable = "on";
            elseif value == 0
                app.CVType.Enable = "off";
                app.Folds.Enable = "off";
                app.Repeats.Enable = "off";                
            end
        end

        % Callback function
        function PredictorModalityButtonGroupSelectionChanged(app, event)
      
        end

        % Callback function
        function PredictorTypeButtonGroupButtonDown(app, event)
            
        end

        % Button pushed function: LoadTextFileButton
        function LoadTextFileButtonPushed(app, event)
            [txt_file, txt_path] = uigetfile({'*.txt' ; '*.'},...
                'Select CSV file containing behavioral data');
            if txt_file == 0
                error('No file selected; please try again');
            else
                app.FileSelText.Value = fullfile(txt_path, txt_file);
            end      
        end

        % Button pushed function: SelectOutputDirectoryButton
        function SelectOutputDirectoryButtonPushed2(app, event)
            [out_dir] = uigetdir(pwd, 'Select directory to save results');
            if out_dir == 0
                error('No directory selected; please try again');
            else
                app.OutDirText.Value = out_dir;
            end           
        end

        % Value changed function: RunPerm
        function RunPermValueChanged(app, event)
            value = app.RunPerm.Value;
            if value == 1
                app.PermIter.Enable = 1;
            else
                app.PermIter.Enable= 0;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 696 453];
            app.UIFigure.Name = 'MATLAB App';

            % Create SelectModelLabel
            app.SelectModelLabel = uilabel(app.UIFigure);
            app.SelectModelLabel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.SelectModelLabel.HorizontalAlignment = 'right';
            app.SelectModelLabel.FontSize = 14;
            app.SelectModelLabel.Position = [9 287 90 22];
            app.SelectModelLabel.Text = 'Select Model:';

            % Create SelectModelDropDown
            app.SelectModelDropDown = uidropdown(app.UIFigure);
            app.SelectModelDropDown.Items = {'plsr', 'pls_da', 'ridge', 'lasso', 'logistic_ridge', 'logistic_lasso', 'linsvr', 'linsvc', 'kernsvr', 'kernsvc', 'rensemble', 'censemble', 'mean', 'mean_score'};
            app.SelectModelDropDown.ValueChangedFcn = createCallbackFcn(app, @SelectModelDropDownValueChanged, true);
            app.SelectModelDropDown.FontSize = 14;
            app.SelectModelDropDown.BackgroundColor = [0.651 0.651 0.651];
            app.SelectModelDropDown.Position = [103 283 100 30];
            app.SelectModelDropDown.Value = 'plsr';

            % Create UseParallelProcessingRecommendedCheckBox
            app.UseParallelProcessingRecommendedCheckBox = uicheckbox(app.UIFigure);
            app.UseParallelProcessingRecommendedCheckBox.Text = 'Use Parallel Processing (Recommended)';
            app.UseParallelProcessingRecommendedCheckBox.FontSize = 14;
            app.UseParallelProcessingRecommendedCheckBox.Position = [215 285 282 22];
            app.UseParallelProcessingRecommendedCheckBox.Value = true;

            % Create StackedModelingAnalysisConfigurationLabel
            app.StackedModelingAnalysisConfigurationLabel = uilabel(app.UIFigure);
            app.StackedModelingAnalysisConfigurationLabel.FontSize = 16;
            app.StackedModelingAnalysisConfigurationLabel.FontWeight = 'bold';
            app.StackedModelingAnalysisConfigurationLabel.Position = [11 427 503 22];
            app.StackedModelingAnalysisConfigurationLabel.Text = 'Stacked Modeling Analysis Configuration';

            % Create RunAnalysisButton
            app.RunAnalysisButton = uibutton(app.UIFigure, 'push');
            app.RunAnalysisButton.ButtonPushedFcn = createCallbackFcn(app, @RunAnalysisButtonPushed, true);
            app.RunAnalysisButton.BackgroundColor = [0.651 0.651 0.651];
            app.RunAnalysisButton.FontSize = 18;
            app.RunAnalysisButton.FontWeight = 'bold';
            app.RunAnalysisButton.Position = [462 47 224 48];
            app.RunAnalysisButton.Text = 'Run Analysis';

            % Create LoadTextFileButton
            app.LoadTextFileButton = uibutton(app.UIFigure, 'push');
            app.LoadTextFileButton.ButtonPushedFcn = createCallbackFcn(app, @LoadTextFileButtonPushed, true);
            app.LoadTextFileButton.BackgroundColor = [0.651 0.651 0.651];
            app.LoadTextFileButton.Position = [7 386 135 20];
            app.LoadTextFileButton.Text = 'Load Text File';

            % Create SelectOutputDirectoryButton
            app.SelectOutputDirectoryButton = uibutton(app.UIFigure, 'push');
            app.SelectOutputDirectoryButton.ButtonPushedFcn = createCallbackFcn(app, @SelectOutputDirectoryButtonPushed2, true);
            app.SelectOutputDirectoryButton.BackgroundColor = [0.651 0.651 0.651];
            app.SelectOutputDirectoryButton.Position = [7 337 135 24];
            app.SelectOutputDirectoryButton.Text = 'Select Output Directory';

            % Create OutDirText
            app.OutDirText = uitextarea(app.UIFigure);
            app.OutDirText.Placeholder = 'No directory selected';
            app.OutDirText.Position = [152 338 533 23];

            % Create FileSelText
            app.FileSelText = uitextarea(app.UIFigure);
            app.FileSelText.Placeholder = 'No file selected';
            app.FileSelText.Position = [152 385 533 23];

            % Create CrossValidationOptionsPanel
            app.CrossValidationOptionsPanel = uipanel(app.UIFigure);
            app.CrossValidationOptionsPanel.Title = 'Cross-Validation Options';
            app.CrossValidationOptionsPanel.FontWeight = 'bold';
            app.CrossValidationOptionsPanel.FontSize = 14;
            app.CrossValidationOptionsPanel.Position = [275 48 180 207];

            % Create PermutationTestingPanel
            app.PermutationTestingPanel = uipanel(app.CrossValidationOptionsPanel);
            app.PermutationTestingPanel.BorderType = 'none';
            app.PermutationTestingPanel.Title = 'Permutation Testing';
            app.PermutationTestingPanel.FontWeight = 'bold';
            app.PermutationTestingPanel.Position = [4 122 143 59];

            % Create RunPerm
            app.RunPerm = uicheckbox(app.PermutationTestingPanel);
            app.RunPerm.ValueChangedFcn = createCallbackFcn(app, @RunPermValueChanged, true);
            app.RunPerm.Text = 'Run Permutations';
            app.RunPerm.Position = [9 9 118 22];

            % Create PermutationsEditField_2Label
            app.PermutationsEditField_2Label = uilabel(app.CrossValidationOptionsPanel);
            app.PermutationsEditField_2Label.HorizontalAlignment = 'right';
            app.PermutationsEditField_2Label.Enable = 'off';
            app.PermutationsEditField_2Label.Position = [6 101 76 22];
            app.PermutationsEditField_2Label.Text = 'Permutations';

            % Create PermIter
            app.PermIter = uieditfield(app.CrossValidationOptionsPanel, 'numeric');
            app.PermIter.Limits = [1 Inf];
            app.PermIter.Enable = 'off';
            app.PermIter.Position = [99 101 37 22];
            app.PermIter.Value = 50;

            % Create Group1EditField_2Label
            app.Group1EditField_2Label = uilabel(app.CrossValidationOptionsPanel);
            app.Group1EditField_2Label.HorizontalAlignment = 'right';
            app.Group1EditField_2Label.Position = [5 40 55 22];
            app.Group1EditField_2Label.Text = 'Group +1';

            % Create PosGroupCost
            app.PosGroupCost = uieditfield(app.CrossValidationOptionsPanel, 'numeric');
            app.PosGroupCost.Limits = [0 Inf];
            app.PosGroupCost.Position = [75 40 25 22];
            app.PosGroupCost.Value = 1;

            % Create Group1EditFieldLabel
            app.Group1EditFieldLabel = uilabel(app.CrossValidationOptionsPanel);
            app.Group1EditFieldLabel.HorizontalAlignment = 'right';
            app.Group1EditFieldLabel.Position = [8 9 52 22];
            app.Group1EditFieldLabel.Text = 'Group -1';

            % Create NegGroupCost
            app.NegGroupCost = uieditfield(app.CrossValidationOptionsPanel, 'numeric');
            app.NegGroupCost.Limits = [0 Inf];
            app.NegGroupCost.Position = [75 9 25 22];
            app.NegGroupCost.Value = 1;

            % Create MisclassificationCostsLabel
            app.MisclassificationCostsLabel = uilabel(app.CrossValidationOptionsPanel);
            app.MisclassificationCostsLabel.FontWeight = 'bold';
            app.MisclassificationCostsLabel.Position = [5 65 138 22];
            app.MisclassificationCostsLabel.Text = 'Misclassification Costs';

            % Create GeneralModelingOptionsPanel
            app.GeneralModelingOptionsPanel = uipanel(app.UIFigure);
            app.GeneralModelingOptionsPanel.Title = 'General Modeling Options';
            app.GeneralModelingOptionsPanel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.GeneralModelingOptionsPanel.FontWeight = 'bold';
            app.GeneralModelingOptionsPanel.FontSize = 14;
            app.GeneralModelingOptionsPanel.Position = [27 48 230 207];

            % Create StandardizeYCheckBox
            app.StandardizeYCheckBox = uicheckbox(app.GeneralModelingOptionsPanel);
            app.StandardizeYCheckBox.Text = 'Standardize Y';
            app.StandardizeYCheckBox.Position = [6 152 97 22];

            % Create StandardizeXCheckBox
            app.StandardizeXCheckBox = uicheckbox(app.GeneralModelingOptionsPanel);
            app.StandardizeXCheckBox.Text = 'Standardize X';
            app.StandardizeXCheckBox.Position = [108 152 97 22];

            % Create HyperParameterTuningPanel
            app.HyperParameterTuningPanel = uipanel(app.GeneralModelingOptionsPanel);
            app.HyperParameterTuningPanel.BorderType = 'none';
            app.HyperParameterTuningPanel.Title = 'Hyper-Parameter Tuning';
            app.HyperParameterTuningPanel.FontWeight = 'bold';
            app.HyperParameterTuningPanel.Position = [6 11 205 134];

            % Create CVStrategyLabel
            app.CVStrategyLabel = uilabel(app.HyperParameterTuningPanel);
            app.CVStrategyLabel.Position = [11 65 70 22];
            app.CVStrategyLabel.Text = 'CV Strategy';

            % Create NFoldsLabel
            app.NFoldsLabel = uilabel(app.HyperParameterTuningPanel);
            app.NFoldsLabel.Position = [11 36 94 22];
            app.NFoldsLabel.Text = 'Number of Folds';

            % Create NRepeatsLabel
            app.NRepeatsLabel = uilabel(app.HyperParameterTuningPanel);
            app.NRepeatsLabel.Position = [11 9 109 22];
            app.NRepeatsLabel.Text = 'Number of Repeats';

            % Create CVType
            app.CVType = uidropdown(app.HyperParameterTuningPanel);
            app.CVType.Items = {'KFold', 'LOOCV'};
            app.CVType.Position = [131 65 69 22];
            app.CVType.Value = 'KFold';

            % Create Folds
            app.Folds = uieditfield(app.HyperParameterTuningPanel, 'numeric');
            app.Folds.Limits = [1 Inf];
            app.Folds.Position = [160 36 40 22];
            app.Folds.Value = 5;

            % Create Repeats
            app.Repeats = uieditfield(app.HyperParameterTuningPanel, 'numeric');
            app.Repeats.Limits = [1 Inf];
            app.Repeats.Position = [160 8 40 22];
            app.Repeats.Value = 5;

            % Create OptimizeHyperparametersCheckBox
            app.OptimizeHyperparametersCheckBox = uicheckbox(app.HyperParameterTuningPanel);
            app.OptimizeHyperparametersCheckBox.Text = 'Optimize Hyper-parameters';
            app.OptimizeHyperparametersCheckBox.Position = [11 89 169 22];
            app.OptimizeHyperparametersCheckBox.Value = true;

            % Create PlottingOptionsPanel
            app.PlottingOptionsPanel = uipanel(app.UIFigure);
            app.PlottingOptionsPanel.Title = 'Plotting Options';
            app.PlottingOptionsPanel.FontWeight = 'bold';
            app.PlottingOptionsPanel.FontSize = 14;
            app.PlottingOptionsPanel.Position = [473 168 191 87];

            % Create ComparisonPlotsPanel
            app.ComparisonPlotsPanel = uipanel(app.PlottingOptionsPanel);
            app.ComparisonPlotsPanel.BorderType = 'none';
            app.ComparisonPlotsPanel.Title = 'Comparison Plots';
            app.ComparisonPlotsPanel.FontWeight = 'bold';
            app.ComparisonPlotsPanel.Position = [4 10 176 51];

            % Create MakePlots
            app.MakePlots = uicheckbox(app.ComparisonPlotsPanel);
            app.MakePlots.Text = 'Make Comparison Plots';
            app.MakePlots.Position = [9 1 168 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = run_stacking_gui

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

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