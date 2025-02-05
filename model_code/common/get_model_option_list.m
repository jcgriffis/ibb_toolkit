function get_model_option_list()

disp('------ Multivariate regression models: ------');
disp('plsr: PLS regression');
disp('ridge: Ridge regression');
disp('lasso: LASSO regression');
disp('rlinsvr: Ridge-regularized linear SVR');
disp('linsvr: Linear SVR');
disp('kernsvr: RBF SVR');
disp('rensemble: Random forest');

disp(' ');
disp('----- Multivariate classification models: ------');
disp('pls_da: PLS linear classifier');
disp('logistic_ridge: Ridge linear classifier');
disp('logistic_lasso: LASSO linear classifier');
disp('rlinsvc: Ride-regularized linear SVC');
disp('linsvc: Linear SVC');
disp('kernsvc: RBF SVC');
disp('censemble: Random forest')

disp(' ');
disp('------ Mass univariate models: ------');
disp('municorr: Peasron correlation');
disp('ttest: T-test (predictors define groups)');
disp('bmunz: Brunner-Munzel (predictors define groups)');
disp('munilr: Logistic regression (outcome defines groups)');
disp('muniolsr: Linear regression')
disp('munimnr: Ordinal regression (ordinal scale outcome)')
disp('prop_sub: Proportional subtraction (outcome defines groups)')

end


