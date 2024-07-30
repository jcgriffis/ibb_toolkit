function [stacked_results, base_results] = create_stacked_models(model_paths, cfg)

% This function takes existing cross-validation predictions from different models, 
% and re-runs the cross-validation to evaluate a stacked model combining the predictions 
% across models into a single model

% Go to model paths and load base model results
base_results.all.pred_y = cell(size(model_paths));
base_results.all.explained = cell(size(model_paths));
base_results.all.r2_ss = cell(size(model_paths));
base_results.pred_y = cell(size(model_paths));
base_results.explained = cell(size(model_paths));
base_results.r2_ss = cell(size(model_paths));
for i = 1:length(model_paths)

    % Load results
    load(fullfile(model_paths{i}, 'cv_results.mat'));

    % Get performance measures collapsed across test folds
    base_results.all.explained{i} = cv_results.all.explained;
    base_results.all.r2_ss{i} = cv_results.all.r2_ss;
    base_results.all.pred_y{i} = cv_results.all.pred_y;    

    % Get fold-level performance measures 
    base_results.explained{i} = cv_results.explained;
    base_results.r2_ss{i} = cv_results.r2_ss;
    base_results.pred_y{i} = cv_results.pred_y;

    % Get relevant data    
    if i == 1
        train_set = cv_results.train_set;
        test_set = cv_results.test_set;
        obs_y = cv_results.all.obs_y(:,1);
        stacked_results.train_set = train_set;
        stacked_results.test_set = test_set;
        stacked_results.obs_y = obs_y; 
    else
        assert(isequal(train_set, cv_results.train_set), 'Error: Training datasets must be the same to allow for model stacking.')
    end

end 

% Set parallel options
if cfg.parallel == 1
    options = statset('UseParallel', true); % parallel processing
else
    options = statset('UseParallel', false); % no parallel processing
end

% Outer loop repeats
stacked_results.explained = zeros(size(train_set, 2), size(train_set, 3));
stacked_results.mse = zeros(size(train_set, 2), size(train_set, 3));
stacked_results.r2_ss = zeros(size(train_set, 2), size(train_set, 3));
stacked_results.corr = zeros(size(train_set, 2), size(train_set, 3));
stacked_results.pred_y = zeros(size(train_set, 1), size(train_set, 2), size(train_set,3)).*NaN;
stacked_results.coeff = zeros(length(model_paths), size(train_set, 2), size(train_set, 3));
stacked_results.coeff_0 = zeros(1, size(train_set, 2), size(train_set, 3));

for i = 1:size(train_set, 2)
    
    % Stack predictors across models
    X_stack = zeros(size(base_results.pred_y{1}, 1), length(base_results.pred_y));
    for k = 1:length(base_results.pred_y)
        % Predictors are predictions from each model for this iteration of repeats and folds
        X_stack(:,k) = base_results.all.pred_y{k}(:,i);          
    end

    % Outer loop folds
    disp(['Running stacked models. Outer loop repeat:' num2str(i) '/' num2str(size(train_set,2))]);
    
    for j = 1:size(train_set, 3)

        % Current training set
        my_train = train_set(:,i,j);
        
        % Current test set
        my_test = test_set(:,i,j);
    
        % Get training and test data
        x_train = X_stack(my_train==1,:);
        x_test = X_stack(my_test==1,:);
        y_train = obs_y(my_train==1);
        y_test = obs_y(my_test==1);

        % Make stacked model
        if strcmp(cfg.model_spec, 'lasso')

            % Train lasso model and optimize hyper-parameters
            [B, fitinfo] = lasso(x_train, y_train, 'Alpha', 1, 'CV', 10, 'MCReps', 10, 'Options', options);
            mdl_idx = fitinfo.IndexMinMSE;
            coef = B(:,mdl_idx);
            coef0 = fitinfo.Intercept(mdl_idx);  
            
            % Get either model predictions or averages of selected features
            if strcmp(cfg.ens_type, 'model_pred')
                y_pred = x_test*coef + coef0;
            elseif strcmp(cfg.ens_type, 'avg_pred')
                y_pred = nanmean(x_test(:, coef~=0),2);
            end
            
            % Get model coefficients
            stacked_results.coeff(:,i,j) = coef;
            stacked_results.coeff_0(:,i,j) = coef0;

        elseif strcmp(cfg.model_spec, 'olsr')
            
            % Fit OLSR model and get predictions
            mdl = fitlm(x_train, y_train);
            y_pred = predict(mdl, x_test);

        elseif strcmp(cfg.model_spec, 'none') 
            y_pred = nanmean(x_test, 2);
        end

        % Get model performance measures for regression
        stacked_results.explained(i,j) = get_explained_variance(y_test, y_pred);
        stacked_results.mse(i,j) = get_mse(y_test, y_pred);
        stacked_results.corr(i,j) = corr(y_test, y_pred);
        stacked_results.r2_ss(i,j) = get_model_r2(y_test, y_pred);
        stacked_results.pred_y(my_test==1,i,j) = y_pred;
    end

end
   
% Get mean performance (across folds)
stacked_results.all.pred_y = nanmean(stacked_results.pred_y,3);
stacked_results.all.explained = zeros(size(stacked_results.all.pred_y,2),1);
stacked_results.all.r2_ss = zeros(size(stacked_results.all.pred_y,2),1);
for i = 1:size(stacked_results.all.pred_y ,2)
    stacked_results.all.explained(i,1) = get_explained_variance(obs_y, stacked_results.all.pred_y(:,i));
    stacked_results.all.r2_ss(i,1) = get_model_r2(obs_y, stacked_results.all.pred_y(:,i));
    stacked_results.all.mse(i,1) = get_mse(obs_y, stacked_results.all.pred_y(:,i));
    stacked_results.all.corr(i,1) = corr(obs_y, stacked_results.all.pred_y(:,i));    
end

end