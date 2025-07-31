function [X, b] = regress_confounds_from_X_cv(X, confounds)

b = zeros(size(X,2),size(confounds,2)+1)';
for i = 1:size(X,2)
    b(:,i) = regress(X(:,i), [ones(size(X,1),1), zscore(confounds)]);
    x_pred = [ones(size(X,1), 1) zscore(confounds)]*b(:,i); % get fitted Y

    % Get Y residuals
    X(:,i) = X(:,i)-x_pred; 
end

end