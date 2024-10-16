function [X_trimmed, X_trim_ind] = trim_X_data(X, bin_thresh, f_thresh)

% Removes predictor variables with < f_thresh observations equal to or greater than bin_thresh
% Joseph Griffis 2023

X_bin = abs(X) >= bin_thresh; % binarize data
X_trim_ind = find(sum(X_bin,1) > round(f_thresh)); % find columns that satisfy required observations 
X_trimmed = X(:,X_trim_ind); % remove columns with too few observations from X

end
