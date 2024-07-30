function model_results = run_perm_mass_uni_lr(X, Y, cfg, model_results)
    
% Permutation tests for mass univariate logistic regressions - takes forever
% Joseph Griffis 2024

% Get permuted Y
n_perm = cfg.perm.n_perm;
perm_Y = zeros(length(Y), n_perm);        
for i = 1:n_perm
    perm_Y(:,i) = Y(randperm(length(Y))');
end    
% glmfit requires Y to be in range [0,1]
if min(perm_Y) == -1 
    perm_Y(perm_Y < 0) = 0;
end
% Preallocate permutation correlation matrix
t_perm = zeros(size(X,2), n_perm);

% Run permutations
if cfg.parallel == 1
    for j = 1:n_perm
        perm_Yj = perm_Y(:,j);
        disp(['Running permutation analysis:' num2str(j)])
        parfor i = 1:length(t_perm)
            [~,~,stats] = glmfit(X(:,i), perm_Yj, 'binomial');
            t_perm(i,j) = stats.t(2);
        end
    end
else
    for j = 1:n_perm
        disp(['Running permutation analysis:' num2str(j)])
        for i = 1:length(t_perm)
            [~,~,stats] = glmfit(X(:,i), perm_Y(:,j), 'binomial');
            t_perm(i,j) = stats.t(2);
        end
    end
end

% Get cFWE thresholds for different numbers of voxels
if cfg.perm.coeff_cfwe == 1
    model_results.perm.cfwe.crit_val_1 = get_wholebrain_cfwe_thresh(t_perm, cfg.fwe_thresh, 1, n_perm);
    model_results.perm.cfwe.crit_val_10 = get_wholebrain_cfwe_thresh(t_perm, cfg.fwe_thresh, 10, n_perm);
    model_results.perm.cfwe.crit_val_50 = get_wholebrain_cfwe_thresh(t_perm, cfg.fwe_thresh, 50, n_perm);
    model_results.perm.cfwe.crit_val_100 = get_wholebrain_cfwe_thresh(t_perm, cfg.fwe_thresh, 100, n_perm);
    model_results.perm.cfwe.crit_val_500 = get_wholebrain_cfwe_thresh(t_perm, cfg.fwe_thresh, 500, n_perm);
    model_results.perm.cfwe.crit_val_1000 = get_wholebrain_cfwe_thresh(t_perm, cfg.fwe_thresh, 1000, n_perm);
end

% Save permutation results if indicated
if cfg.perm.save_perm_results == 1
    if isfile('perm_results.mat')
       save('perm_results.mat', 't_perm', '-append');
    else
       save('perm_results.mat', 't_perm', '-v7.3');
    end
end
clear perm_results

end
