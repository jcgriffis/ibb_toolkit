function X = regress_confounds_from_X(X, confounds)

for i = 1:size(X,2)
    b = regress(X(:,i), [ones(size(X,1),1), zscore(confounds)]);
    yhat = [ones(size(X,1), 1) zscore(confounds)]*b; % get fitted Y

    % Get Y residuals
    y_resid = X(:,i)-yhat;
    X(:,i) = y_resid;
end

end