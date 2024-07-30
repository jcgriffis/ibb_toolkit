function [y_resid, r2, b] = regress_confounds_from_Y(Y, confounds)

% Perform a linear regression to regress out mean effects of confound variables from response; residuals are defined as new response variable
% Joseph Griffis 2023

% Fit a linear regression predicting Y based on lesion volume
b = regress(Y, [ones(length(Y), 1), confounds]);
yhat = [ones(length(Y), 1) confounds]*b; % get fitted Y

% Get Y residuals
y_resid = Y-yhat;

% Get variance explained by confound model
r2 = get_explained_variance(Y, yhat);

end