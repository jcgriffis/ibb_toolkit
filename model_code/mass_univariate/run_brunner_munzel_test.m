function model_results = run_brunner_munzel_test(X, Y, cfg, model_results)

% Run Brunner-Munzel test
% Joseph Griffis 2024

% Test direction
if strcmp(cfg.direction, 'pos')
    sidedness = 3;
elseif strcmp(cfg.direction, 'neg')
    sidedness = 2;
else
    sidedness = 1;
end

% Preallocate values
pval = zeros(size(X,2),1);
bmstat = zeros(size(X,2),1);

% Run test
parfor i = 1:size(X,2)
    bmstats = approxbrunnermunzel(Y(X(:,i)==0), Y(X(:,i)==1), sidedness);
    bmstat(i) = bmstats.stat;
    pval(i) = bmstats.pval;
end

% Get voxel level multiple comparisons 
model_results.coeff_vfdr_pvals = mafdr(pval, 'BHFDR', true);
model_results.coeff_vfwe_pvals = bonf_holm(pval, cfg.fwe_thresh);
model_results.coeff = bmstat;

% Run permutations
if cfg.permutation == 1

    % number of permutation iterations
    n_perm = cfg.perm.n_perm;
    
    % Preallocate permutation result matrix
    bm_perm = zeros(size(X,2), n_perm);

    if cfg.parallel == 1
        for i = 1:n_perm
            disp(['Running permutation analysis: ' num2str(i)]);
            perm_Y = Y(randperm(length(Y))');
            parfor j = 1:size(X,2)
                permstat = approxbrunnermunzel(perm_Y(X(:,j)==0), perm_Y(X(:,j)==1), sidedness);
                bm_perm(j,i) = permstat.stat; 
            end
        end
    else
        for i = 1:n_perm
            disp(['Running permutation analysis: ' num2str(i)]);
            perm_Y = Y(randperm(length(Y))');
            for j = 1:size(X,2)
                permstat = approxbrunnermunzel(perm_Y(X(:,j)==0), perm_Y(X(:,j)==1), sidedness);
                bm_perm(j,i) = permstat.stat;             
            end
        end
    end
    
    % Compute p-values
    disp('Computing permutation p-values...');
    model_results = get_permutation_results(model_results, bm_perm, cfg);
    
    disp('Finished permutation analysis');

    % Save permutation results if indicated
    if cfg.perm.save_perm_results == 1
        if isfile('perm_results.mat')
           save('perm_results.mat', 'bm_perm', '-append');
        else
           save('perm_results.mat', 'bm_perm', '-v7.3');
        end
    end
    clear perm_results
end

% Store Y
model_results.obs_y = Y;

end

