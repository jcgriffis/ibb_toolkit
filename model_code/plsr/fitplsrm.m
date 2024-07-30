function [XS, betas, yhat, r2, vip_score, mse] = fitplsrm(X, Y, opt_k)

% Fits a PLSR rgression model with specified dimensionality and computes R-squared
% Joseph Griffis 2023
    
% Fit PLS regression
[XL,yl,XS,~,betas,~,~,stats] = plsregress(X, Y, opt_k);

% Get vip scores
W0 = stats.W ./ sqrt(sum(stats.W.^2,1));
p = size(XL,1);
sumSq = sum(XS.^2,1).*sum(yl.^2,1);
vip_score = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));

% Get predictions and R-squared
yhat = [ones(length(Y), 1) X]*betas; % get fitted Y
r2 = get_model_r2(Y, yhat); % r-squared
mse = get_mse(Y, yhat);

end