function [r2] = get_explained_variance(obs_y, pred_y)

% Compute R-squared using explained variance method
% Joseph Griffis 2024
numerator = var(obs_y - pred_y);
denominator = var(obs_y);
r2 = 1 - (numerator/denominator);
if isinf(r2)
    r2 = NaN;
end

end