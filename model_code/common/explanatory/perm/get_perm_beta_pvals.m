function p_vals = get_perm_beta_pvals(b_real, b_perm, n_perm)

% Compute PLSR beta p-values based on permutation results
% Joseph Griffis 2023

% Make sure permutation betas are oriented correctly
if size(b_perm,2) ~= n_perm
    b_perm = b_perm';
end

%%%%  Find negative differences between real and permuted absolute beta magnitudes
% i.e., if a permuted beta is a larger weight than a true beta, then the difference will be negative
% We can count the number of negative differences to generate a p-value

% Get differences in absolute beta magnitudes
p_vals = (sum((abs(b_real) - abs(b_perm)) < 0, 2)+1)./(n_perm+1);

end
