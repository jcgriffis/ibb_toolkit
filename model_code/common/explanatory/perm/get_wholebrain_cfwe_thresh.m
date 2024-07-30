function crit_val = get_wholebrain_cfwe_thresh(beta_perm, alpha, n_vox, n_perm)

% Apply whole-brain cFWE thresholding based on max permutation statistics (Mirman et al., 2018 - Neuropsychologia)
% Joseph Griffis 2024

% Make sure permutation betas are oriented correctly
if size(beta_perm,2) ~= n_perm
    beta_perm = beta_perm';
end

% Sort permutation betas
beta_perm = sort(abs(beta_perm), 'descend');

if n_vox > size(beta_perm,1)
    crit_val = NaN;
else
    % Get threshold for selected alpha
    crit_val = prctile(beta_perm(n_vox,:), 100.*(1-alpha));
end

end