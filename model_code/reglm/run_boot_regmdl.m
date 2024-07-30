function [out] = run_boot_regmdl(X, Y, lambda, learner, reg_type, standardize)

% Run bootstrap analyses
% Joseph Griffis 2024

% Get model params
if learner == 1
    learner = 'svm';
else
    learner = 'leastsquares';
end
if reg_type == 1
    reg_type = 'lasso';
else
    reg_type = 'ridge';
end

if standardize == 1
    [X] = normalize(X);
    X(isnan(X))=0; % Since normalize will cause columns to become NaN if they are all 0
elseif standardize == 2
    Y = normalize(Y);
elseif standardize == 3
    [X] = normalize(X);
    X(isnan(X))=0; % Since normalize will cause columns to become NaN if they are all 0   
    Y = normalize(Y);
end

% Fit final model
mdl = fitrlinear(X,Y,...
    'ObservationsIn', 'rows',...
    'Regularization', reg_type,...
    'Lambda', lambda,...
    'Learner', learner,...
    'BetaTolerance',0);

% Get betas
betas = mdl.Beta;

% Get final model R-squared
y_pred = predict(mdl, X, 'ObservationsIn', 'rows');
r2 = get_model_r2(Y, y_pred);
out = [r2; betas];

end