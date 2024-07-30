function [model_results] = run_perm_svc(X, Y, cfg, model_results)

% Use K-fold cross-validation to tune regularized regression hyper-parameters
% and fit final model to full dataset using permutation testing 
% Joseph Griffis 2023

% Permutation variables
n_perm = cfg.perm.n_perm;
roc_auc = zeros(n_perm, 1);
if cfg.perm.coeff_p == 1 || cfg.perm.coeff_cfwe == 1
    beta_perm = zeros(size(X,2), n_perm);
    perm_coeff = 1;
else
    perm_coeff = 0;
end

% Avoid dot indexing in parallel processing
C = model_results.C;
gamma = model_results.gamma;
kernel = cfg.kernel;
S.ClassNames = [1, -1];
S.ClassificationCosts = cfg.cost;

% Fit null models
if cfg.parallel == 1

    parfor i = 1:n_perm           

        % Fit final model with cross-validation optimized hyper-parameters
        disp(['Permutation ' num2str(i)]);
        perm_Y = Y(randperm(length(Y)));  
        
        alpha = [];
        while isempty(alpha)
            null_model = fitcsvm(X,perm_Y,'KernelFunction', kernel,...
                'BoxConstraint', C,...
                'KernelScale', gamma, ...
                'Cost', S);
            if isnan(null_model.Bias)
                alpha = [];
            else
                alpha = null_model.Alpha;
            end
        end

        % Get predictions
        [~, pred_score] = predict(null_model, X);

        % Compute area under ROC curve
        roc_obj = rocmetrics(perm_Y, pred_score(:,2), 1);
        roc_auc(i) = roc_obj.AUC;
        if perm_coeff == 1
            if strcmp(kernel, 'rbf')
                beta_perm(:,i) = (null_model.Alpha.'*(null_model.SupportVectors.*null_model.SupportVectorLabels)); 
            else
                beta_perm(:,i) = null_model.Beta;
            end
        end
    end
else
    for i = 1:n_perm
        
        % Fit final model with cross-validation optimized hyper-parameters
        disp(['Permutation ' num2str(i)]);
        perm_Y = Y(randperm(length(Y)));  
        
        alpha = [];
        while isempty(alpha)
            null_model = fitcsvm(X,perm_Y,'KernelFunction', kernel,...
                'BoxConstraint', C,...
                'KernelScale', gamma, ...
                'Cost', S);
            if isnan(null_model.Bias)
                alpha = [];
            else
                alpha = null_model.Alpha;
            end
        end

        % Get predictions
        [~, pred_score] = predict(null_model, X);

        % Compute area under ROC curve
        roc_obj = rocmetrics(perm_Y, pred_score(:,2), 1);
        roc_auc(i) = roc_obj.AUC;
        if perm_coeff == 1
            if strcmp(kernel, 'rbf')
                beta_perm(:,i) = (null_model.Alpha.'*(null_model.SupportVectors.*null_model.SupportVectorLabels)); 
            else
                beta_perm(:,i) = null_model.Beta;
            end
        end
        clear null_model
    end
end

% Get p-values
model_results.perm.model_pval = (numel(find(roc_auc >= model_results.roc_auc))+1)./(n_perm+1);
if perm_coeff == 1
    model_results = get_permutation_results(model_results, beta_perm, cfg);
end
disp('Finished permutation analysis');

% Save permutation results if indicated
if cfg.perm.save_perm_results == 1
    if isfile('perm_results.mat')
       if exist('beta_perm', 'var')
           save('perm_results.mat', 'roc_auc',  'beta_perm', '-append');
       else
           save('perm_results.mat', 'roc_auc', '-append');
       end           
    else
        if exist('beta_perm', 'var')
           save('perm_results.mat', 'roc_auc',  'beta_perm', '-v7.3');
        else
           save('perm_results.mat', 'roc_auc', '-v7.3');
        end
    end
end

end