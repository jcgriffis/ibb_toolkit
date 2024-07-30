function model_results = get_permutation_coefficient_pvalues(model_results, field, perm_stat, fwe_thresh, n_perm)

% Apply whole-brain cFWE thresholding for v=[1,10,50,100,500,1000] (Mirman et al., 2018 - Neuropsychologia)
% Joseph Griffis 2024

disp('Permutation tests complete, computing whole-brain continuous FWE p-value thresholds...');
model_results.perm.(field).crit_val_1 = get_wholebrain_cfwe_thresh(perm_stat, fwe_thresh, 1, n_perm);    
model_results.perm.(field).crit_val_10 = get_wholebrain_cfwe_thresh(perm_stat, fwe_thresh, 10, n_perm);
model_results.perm.(field).crit_val_50 = get_wholebrain_cfwe_thresh(perm_stat, fwe_thresh, 50, n_perm);
model_results.perm.(field).crit_val_100 = get_wholebrain_cfwe_thresh(perm_stat, fwe_thresh, 100, n_perm);
model_results.perm.(field).crit_val_500 = get_wholebrain_cfwe_thresh(perm_stat, fwe_thresh, 500, n_perm);
model_results.perm.(field).crit_val_1000 = get_wholebrain_cfwe_thresh(perm_stat, fwe_thresh, 1000, n_perm);

end
