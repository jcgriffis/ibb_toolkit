function [model_results] = run_perm_svr(X, Y, cfg, model_results)

% Use K-fold cross-validation to tune regularized regression hyper-parameters
% and fit final model to full dataset using permutation testing 
% Joseph Griffis 2023

% Permutation variables
n_perm = cfg.perm.n_perm;
mse_perm = zeros(1, n_perm);
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

% Standardization
if isfield(model_results, 'Sy')
    standardize_y = 1;
    Sy = model_results.Sy;
    Cy = model_results.Cy;
else 
    standardize_y = 0;
    Sy = [];
    Cy = [];
end

% Fit null models
if cfg.parallel == 1

    parfor i = 1:n_perm        
        % Fit final model with cross-validation optimized hyper-parameters
        disp(['Permutation ' num2str(i)]);
        perm_Y = Y(randperm(length(Y)));       
        alpha = [];
        while isempty(alpha)
            null_model = fitrsvm(X,perm_Y,'KernelFunction', kernel,...
                'BoxConstraint', C,...
                'KernelScale', gamma);
            if isnan(null_model.Bias)
                alpha = [];
            else
                alpha = null_model.Alpha;
            end
        end
        
        % Get predictions
        pred_y = predict(null_model, X);
        if standardize_y == 1
            pred_y = (pred_y .* Sy) + Cy;
            perm_Y = (perm_Y .* Sy) + Cy;
        end                
        mse_perm(i) = get_mse(perm_Y, pred_y);
        if perm_coeff == 1
            if strcmp(kernel, 'rbf')
                beta_perm(:,i) = (null_model.Alpha.'*null_model.SupportVectors); 
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
            null_model = fitrsvm(X,perm_Y,'KernelFunction', kernel,...
                'BoxConstraint', C,...
                'KernelScale', gamma);
            if isnan(null_model.Bias)
                alpha = [];
            else
                alpha = null_model.Alpha;
            end
        end
        
        % Get predictions
        pred_y = predict(null_model, X);
        if standardize_y == 1
            pred_y = (pred_y .* Sy) + Cy;
            perm_Y = (perm_Y .* Sy) + Cy;
        end                
        mse_perm(i) = get_mse(perm_Y, pred_y);
        if perm_coeff == 1
            if strcmp(kernel, 'rbf')
                beta_perm(:,i) = (null_model.Alpha.'*null_model.SupportVectors); 
            else
                beta_perm(:,i) = null_model.Beta;
            end
        end
        clear null_model
    end
end

% Get p-values
model_results.perm.model_pval = (numel(find(round(mse_perm,4) <= round(model_results.mse,4)))+1)./(n_perm+1);
if perm_coeff == 1
    model_results = get_permutation_results(model_results, beta_perm, cfg);
end
disp('Finished permutation analysis');

% Save permutation results if indicated
if cfg.perm.save_perm_results == 1
    if isfile('perm_results.mat')
       if exist('beta_perm', 'var')
           save('perm_results.mat', 'mse_perm',  'beta_perm', '-append');
       else
           save('perm_results.mat', 'mse_perm', '-append');
       end           
    else
        if exist('beta_perm', 'var')
           save('perm_results.mat', 'mse_perm', 'beta_perm', '-v7.3');
        else
           save('perm_results.mat', 'mse_perm', '-v7.3');
        end
    end
end

end