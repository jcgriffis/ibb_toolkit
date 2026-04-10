function model_results = run_perm_mass_uni_mnr(X, Y, cfg, model_results)
    
% Permutation tests for mass univariate ordinal logistic regressions
% The approach is based on Potter et al., 2005 - Statistics in Medicine
% Joseph Griffis 2024

% Check version
my_version = version('-release');
my_version = str2double(my_version(1:end-1));

% Get permuted Y
n_perm = cfg.perm.n_perm;
perm_order = zeros(length(Y), n_perm);
for i = 1:n_perm
    perm_order(:,i) = randperm(length(Y))';
end    

% Preallocate permutation correlation matrix
t_perm = zeros(size(X,2), n_perm);

% Define confounds (parfor doesn't like dot indexing)
confounds = cfg.confounds;
n_levels = numel(unique(Y));

% Run permutations
if cfg.parallel == 1
    for j = 1:n_perm
        disp(['Running permutation analysis:' num2str(j)])
        my_perm = perm_order(:,j);
        confound_des = [ones(size(X,1),1), confounds];
        if my_version >= 2023
            parfor i = 1:size(t_perm,1)
                warning('off', 'MATLAB:nearlySingularMatrix');        
                warning('off', 'stats:mnrfit:IterOrEvalLimit'); 
                [~,~,res] = regress(X(:,i), confound_des);                
                mdl = fitmnr([res(my_perm), confounds], Y, 'model', 'ordinal');
                t_perm(i,j) = mdl.Coefficients.tStat(n_levels);
            end        
        else
            parfor i = 1:size(t_perm,1)        
                warning('off', 'MATLAB:nearlySingularMatrix');        
                warning('off', 'stats:mnrfit:IterOrEvalLimit');        
                [~,~,res] = regress(X(:,i), confound_des);                                
                [~,~,stats] = mnrfit([res(my_perm), confounds], Y, 'model', 'ordinal');
                t_perm(i,j) = stats.t(n_levels);
            end
        end
    end
else
    for j = 1:n_perm
        disp(['Running permutation analysis:' num2str(j)])
        my_perm = perm_order(:,j);
        confound_des = [ones(size(X,1),1), confounds];
        if my_version >= 2023
            for i = 1:size(t_perm,1)
                warning('off', 'MATLAB:nearlySingularMatrix');        
                warning('off', 'stats:mnrfit:IterOrEvalLimit'); 
                [~,~,res] = regress(X(:,i), confound_des);                
                mdl = fitmnr([res(my_perm), confounds], Y, 'model', 'ordinal');
                t_perm(i,j) = mdl.Coefficients.tStat(n_levels);
            end        
        else
            for i = 1:size(t_perm,1)        
                warning('off', 'MATLAB:nearlySingularMatrix');        
                warning('off', 'stats:mnrfit:IterOrEvalLimit');        
                [~,~,res] = regress(X(:,i), confound_des);                                
                [~,~,stats] = mnrfit([res(my_perm), confounds], Y, 'model', 'ordinal');
                t_perm(i,j) = stats.t(n_levels);
            end
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
