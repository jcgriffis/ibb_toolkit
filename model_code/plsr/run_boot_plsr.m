% This function serves as input to the bootci anonymous function.
% Outputs are a matrix of CIs for model fit and betas. 
% Beta CIs don't include the
% intercept.

% Joseph Griffis, 2018, Washington University in St. Louis    

function [out] = run_boot_plsr(X,Y,k, standardize, method, type)

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

    [~,~,~,~,BETA,~,~,~]  = plsregress(X,Y,k); % fit fixed effects model with opt_k   
    
    % Get predictions and R-squared
    yhat = [ones(length(Y), 1) X]*BETA; % get fitted Y
    gof = get_model_r2(Y, yhat);
    
    % Outputs
    betas = BETA(2:end); % get non-intercept betas
    out = [gof;  betas]; % output matrix
end


