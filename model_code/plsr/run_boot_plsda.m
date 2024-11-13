% This function serves as input to the bootci anonymous function.
% Outputs are a matrix of CIs for model fit and betas. 
% Beta CIs don't include the
% intercept.

% Joseph Griffis, 2018, Washington University in St. Louis    

function [out] = run_boot_plsda(X,Y,k, cost, standardize, method, type)

    if standardize == 1 || standardize == 3
        X = normalize(X, method, type);
        X(isnan(X))=0; % Since normalize will cause columns to become NaN if they are all 0
    end

    [~,~,~,~,BETA,~,~,~]  = plsregress(X,Y,k); % fit fixed effects model with opt_k   
    
    % Get predictions and R-squared
    yhat = [ones(length(Y), 1) X]*BETA; % get fitted Y

    % Get ROC AUC 
    roc_obj = rocmetrics(Y, yhat, 1);    
    gof(4,1) = roc_obj.AUC;    
   
    % Get optimal threshold
    opt_thresh = get_optimal_threshold(Y, roc_obj, cost);

    % Convert to labels
    pred_y = zeros(length(yhat),1);
    pred_y(yhat >= opt_thresh) = 1;
    pred_y(yhat < opt_thresh) = -1;
    
    % Compute accuracy for each group and for full sample
    gof(2,1) = numel(intersect(find(pred_y==1), find(Y==1)))./numel(find(Y==1));
    gof(1,1) = numel(intersect(find(pred_y==-1), find(Y==-1)))./numel(find(Y==-1));
    gof(3,1) = numel(find(pred_y==Y))./numel(find(Y));

    % Outputs
    betas = BETA(2:end); % get non-intercept betas
    out = [gof;  betas]; % output matrix
end

