function model_results = run_perm_mass_uni_olsr(X, Y, cfg, model_results)
    
% Permutation tests for mass univariate OLS regressions - takes forever
% Joseph Griffis 2024

% Residualize Y if necessary for Freedman-Lane (Winkler et al., 2014)
if ~isempty(cfg.confounds)
    mdl = fitlm(cfg.confounds, Y);
    Y = mdl.Residuals.Raw;
end

% Get permuted Y
n_perm = cfg.perm.n_perm;
perm_Y = zeros(length(Y), n_perm);        
for i = 1:n_perm
    perm_Y(:,i) = Y(randperm(length(Y))');
end    

% Preallocate permutation matrix
t_perm = zeros(size(X,2), n_perm);

% Get confounds
confounds = cfg.confounds;

% Run permutations
if cfg.parallel == 1
    for j = 1:n_perm
        perm_Yj = perm_Y(:,j);
        disp(['Running permutation analysis:' num2str(j)])
        parfor i = 1:length(t_perm)
            warning('off', 'stats:glmfit:IllConditioned');     
            warning('off', 'stats:glmfit:IterationLimit');            
            [~,~,stats] = glmfit([X(:,i), confounds], perm_Yj, 'normal');
            t_perm(i,j) = stats.t(2);
        end
    end
else
    warning('off', 'stats:glmfit:IllConditioned');     
    warning('off', 'stats:glmfit:IterationLimit');    
    for j = 1:n_perm
        disp(['Running permutation analysis:' num2str(j)])
        for i = 1:length(t_perm)
            [~,~,stats] = glmfit([X(:,i), confounds], perm_Y(:,j), 'normal');
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
