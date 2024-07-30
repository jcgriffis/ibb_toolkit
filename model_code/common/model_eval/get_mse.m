function mse = get_mse(Y, y_hat)

% compute MSE
% Joseph Griffis 2023
sq_err = (Y-y_hat).^2;
mse = mean(sq_err);

end