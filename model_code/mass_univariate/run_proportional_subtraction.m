function model_results = run_proportional_subtraction(X, Y, cfg, model_results)

% Get proportional overlap for each group
model_results.overlap_g0 = sum(X(Y==-1,:),1) ./ numel(Y(Y==-1));
model_results.overlap_g1 = sum(X(Y==1,:),1) ./ numel(Y(Y==1));

% Compute voxel-wise subtraction
model_results.g1_minus_g0 = model_results.overlap_g1 - model_results.overlap_g0;

% Run permutations
if cfg.permutation == 1

    % Get permuted Y
    n_perm = cfg.perm.n_perm;
    
    % Preallocate permutation result matrix
    dif_perm = zeros(size(X,2), n_perm);    

    % Run permutations
    if cfg.parallel == 1
        parfor i = 1:n_perm
            disp(['Running permutation analysis: ' num2str(i)]);
            perm_Y = Y(randperm(length(Y))');
            overlap_g0 = sum(X(perm_Y==-1,:),1) ./ numel(perm_Y(perm_Y==-1));
            overlap_g1 = sum(X(perm_Y==1,:),1) ./ numel(perm_Y(perm_Y==1));
            dif_perm(:,i) = overlap_g1 - overlap_g0;
        end
    else
       for i = 1:n_perm
           disp(['Running permutation analysis: ' num2str(i)]);
           perm_Y = Y(randperm(length(Y))');
           overlap_g0 = sum(X(perm_Y==-1,:),1) ./ numel(perm_Y(perm_Y==-1));
           overlap_g1 = sum(X(perm_Y==1,:),1) ./ numel(perm_Y(perm_Y==1));
           dif_perm(:,i) = overlap_g1 - overlap_g0;
       end     
    end
    
    % Compute p-values
    disp('Computing permutation p-values...');
    model_results = get_permutation_results(model_results, dif_perm, cfg);

    % Save permutation results if indicated
    if cfg.perm.save_perm_results == 1
        if isfile('perm_results.mat')
           save('perm_results.mat', 'dif_perm', '-append');
        else
           save('perm_results.mat', 'dif_perm', '-v7.3');
        end
    end
    clear perm_results
    
    disp('Finished permutation analysis');

end

% Save observed Y
model_results.obs_y = Y;

end