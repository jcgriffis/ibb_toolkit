function [r2] = get_model_r2_train(obs_y, pred_y, train_y)

% Get model R-squared using sum of squares method
% Joseph Griffis 2023 
numerator = sum((obs_y-pred_y).^2); % get sum of squared residuals
denominator = sum((obs_y-mean(train_y)).^2); % get total sum of squares 
r2 = 1-(numerator/denominator); % Get r-squared

end