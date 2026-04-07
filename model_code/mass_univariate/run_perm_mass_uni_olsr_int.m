function model_results = run_perm_mass_uni_olsr_int(X, Y, cfg, model_results)
    
% Permutation tests for mass univariate OLS regressions, uses Freedman-Lane (Winkler et al., 2014)
% Joseph Griffis 2024

% Get permuted Y
n_perm = cfg.perm.n_perm;
perm_order = zeros(length(Y), n_perm);        
for i = 1:n_perm
    perm_order(:,i) = randperm(length(Y))';
end    

% Preallocate permutation matrix
t_perm = zeros(size(X,2), n_perm);

% Get confounds
confounds = cfg.confounds;
int_term = confounds(:,cfg.int_term);

% Run permutations
if cfg.parallel == 1
    for j = 1:n_perm
        my_perm = perm_order(:,j);
        confound_des = [ones(size(X,1),1), confounds];        
        disp(['Running permutation analysis:' num2str(j)])
        parfor i = 1:size(t_perm,1)
            warning('off', 'stats:glmfit:IllConditioned');     
            warning('off', 'stats:glmfit:IterationLimit');       
            res = regress(Y, [confound_des, X(:,i)]);      
            [~,~,stats] = glmfit([X(:,i).*int_term, X(:,i), confounds], res(my_perm), 'normal');
            t_perm(i,j) = stats.t(2);
        end
    end
else
    for j = 1:n_perm
        my_perm = perm_order(:,j);
        confound_des = [ones(size(X,1),1), confounds];        
        disp(['Running permutation analysis:' num2str(j)])
        for i = 1:size(t_perm,1)
            warning('off', 'stats:glmfit:IllConditioned');     
            warning('off', 'stats:glmfit:IterationLimit');       
            res = regress(Y, [confound_des, X(:,i)]);      
            [~,~,stats] = glmfit([X(:,i).*int_term, X(:,i), confounds], res(my_perm), 'normal');
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
