function model_results = run_ttests(X, Y, cfg, model_results)

% Run mass-univariate t-tests - limited functionality
% Joseph Griffis 2024

% Preallocate values
pval = zeros(size(X,2),1);
tstat = zeros(size(X,2),1);

% Run test
vartype = cfg.vartype;
if cfg.parallel == 1
    parfor i = 1:size(X,2)
        [~, pval(i), ~, stats] = ttest2(Y(X(:,i)==1), Y(X(:,i)==0), 'vartype', vartype);
        tstat(i) = stats.tstat;
    end
else
   for i = 1:size(X,2)
        [~, pval(i), ~, stats] = ttest2(Y(X(:,i)==1), Y(X(:,i)==0), 'vartype', vartype);
        tstat(i) = stats.tstat;
   end
end

% Get voxel level multiple comparisons 
model_results.coeff_vfdr_pvals = mafdr(pval, 'BHFDR', true);
model_results.coeff_vfwe_pvals = bonf_holm(pval, cfg.fwe_thresh);
model_results.coeff_pvals = pval;
model_results.tstat = tstat;

% Get permuted Y
n_perm = cfg.perm.n_perm;

% Preallocate permutation correlation matrix
t_perm = zeros(size(X,2), n_perm);

% Run permutations
if cfg.permutation == 1
    if cfg.parallel == 1
        for i = 1:n_perm
            disp(['Running permutation analysis: ' num2str(i)]);
            perm_Y = Y(randperm(length(Y))');
            parfor j = 1:size(X,2)
                [~, ~, ~, stats] = ttest2(perm_Y(X(:,j)==1), perm_Y(X(:,j)==0), 'vartype', 'unequal');
                t_perm(j,i) = stats.tstat; 
            end
        end
    else
        for i = 1:n_perm
            disp(['Running permutation analysis: ' num2str(i)]);
            perm_Y = Y(randperm(length(Y))');
            for j = 1:size(X,2)
                [~, ~, ~, stats] = ttest2(perm_Y(X(:,j)==1), perm_Y(X(:,j)==0), 'vartype', 'unequal');
                t_perm(j,i) = stats.tstat;           
            end
        end
    end
    
    % Compute p-values
    disp('Computing permutation p-values...');
    model_results = get_permutation_results(model_results, t_perm, cfg);

% Save permutation results if indicated
if cfg.perm.save_perm_results == 1
    if isfile('perm_results.mat')
       save('perm_results.mat', 't_perm', '-append');
    else
       save('perm_results.mat', 't_perm', '-v7.3');
    end
end
clear perm_results

disp('Finished permutation analysis');

end

% Save observed Y
model_results.obs_y = Y;

end