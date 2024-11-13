function [out] = run_boot_svr(X, Y, C, gamma, kernel, standardize, method, type)

% Run bootstrap analyses
% Joseph Griffis 2024

% Set kernel function
if kernel == 0
    kfun = 'linear';
else
    kfun = 'rbf';
end

if standardize == 1
    [X] = normalize(X, method, type);
    X(isnan(X))=0; % Since normalize will cause columns to become NaN if they are all 0   
elseif standardize == 2
    Y = normalize(Y, method, type);
elseif standardize == 3
    [X] = normalize(X, method, type);
    X(isnan(X))=0; % Since normalize will cause columns to become NaN if they are all 0   
    Y = normalize(Y, method, type);
end


alpha = [];
while isempty(alpha)
    % Fit final model with cross-validation optimized hyper-parameters
    mdl = fitrsvm(X,Y,'KernelFunction', kfun,...
    'BoxConstraint', C,...
    'KernelScale', gamma);
    if isnan(mdl.Bias)
        alpha = [];
    else
        alpha = mdl.Alpha;
    end
end

% Compute betas using sensitivity mapping method (Zhang et al., 2014 - Human Brain Mapping)
if strcmp(kfun, 'rbf')
    betas(:,1) = (mdl.Alpha.'*mdl.SupportVectors);
else
    betas(:,1) = mdl.Beta;
end

% Get final model R-squared
y_pred = predict(mdl, X);
r2 = get_model_r2(Y, y_pred);
out = [r2; betas];

end