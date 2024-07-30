function [opt_k, y_mse] = run_loocv_plsr(X,Y,cfg)
% This function identifies the optimal number of PLSR components
% Inputs:
% X: predictor matrix (Nxk)
% Y: response vector (Nx1)
% cfg: structure containing parameters/options
% Joseph Griffis 2023

[rows_x, ~] = size(X); % get size of predictor matrix 
    
% Partition data 
c = cvpartition(rows_x, 'LeaveOut');

% Max number of components
if isfield(cfg.hp_opt, 'n_comp')
    n_comp = cfg.hp_opt.n_comp;
else
    n_comp = min(size(X))-1;
end

% Parallel flag
if cfg.parallel == 1
    options = statset('UseParallel',true);
else
    options = statset('UseParallel',false);
end

% Fit PLS model with Kfold cross-validation
[~,~,~,~,~,~, MSE] = plsregress(X, Y, n_comp, 'cv', c, 'Options', options);

% Get MSE for outcome
y_mse = MSE(2,2:end);

% Find last component that decreased CV MSE
opt_k = find(diff(y_mse) > 0,1);

end