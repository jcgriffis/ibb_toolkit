function [model_results] = run_perm_plsda(X, Y, cfg, model_results)

% Run permutation testing for PLSDA
% Joseph Griffis 2024

% preallocate for permutation
n_perm = cfg.perm.n_perm;
roc_auc = zeros(n_perm, 1);
n_perm = cfg.perm.n_perm;
if cfg.perm.coeff_p == 1 || cfg.perm.coeff_cfwe == 1
    beta_perm = zeros(size(X,2), n_perm);
    perm_coeff = 1;
else
    perm_coeff = 0;
end

% Parallel flag
if cfg.parallel == 1
    options = statset('UseParallel',true);
else
    options = statset('UseParallel',false);
end

% Run permutation tests
n_comp = model_results.opt_k; % model dimensionality from CV
disp(['Starting permutation tests with ' num2str(n_perm) ' permutations...']);
if cfg.parallel == 1
    parfor i = 1:n_perm
        disp(['Permutation Test ' num2str(i) '/' num2str(n_perm)]);
        
        % Fit null model
        perm_Y = Y(randperm(length(Y)));
        [~,~,~,~,betas,~,~,~] = plsregress(X, perm_Y, n_comp, 'options', options); % fit fixed effects model with opt_k
        pred_y = [ones(length(perm_Y), 1) X]*betas; % get fitted perm_Y
        
        % Get ROC AUC
        roc_obj = rocmetrics(perm_Y, pred_y, 1);
        roc_auc(i) = roc_obj.AUC;
        
        % Get betas
        if perm_coeff == 1
            beta_perm(:,i) = betas(2:end);
        end
    end
else
    for i = 1:n_perm
        disp(['Permutation Test ' num2str(i) '/' num2str(n_perm)]);
        
        % Fit null model
        perm_Y = Y(randperm(length(Y)));
        [~,~,~,~,betas,~,~,~] = plsregress(X, perm_Y, n_comp, 'options', options); % fit fixed effects model with opt_k
        pred_y = [ones(length(perm_Y), 1) X]*betas; % get fitted perm_Y
        
        % Get ROC AUC
        roc_obj = rocmetrics(perm_Y, pred_y, 1);
        roc_auc(i) = roc_obj.AUC;
               
        % Get betas
        if perm_coeff == 1
            beta_perm(:,i) = betas(2:end);
        end
    end
end

% Get p-values
model_results.perm.model_pval = (numel(find(roc_auc >= model_results.roc_auc))+1)./(n_perm+1);
if perm_coeff == 1
    model_results = get_permutation_results(model_results, beta_perm, cfg);
end

disp('Finished permutation analysis');

% Save permutation results if indicated
if cfg.perm.save_perm_results == 1
    if isfile('perm_results.mat')
       if exist('beta_perm', 'var')
           save('perm_results.mat', 'roc_auc',  'beta_perm', '-append');
       else
           save('perm_results.mat', 'roc_auc', '-append');
       end           
    else
        if exist('beta_perm', 'var')
           save('perm_results.mat', 'roc_auc',  'beta_perm', '-v7.3');
        else
           save('perm_results.mat', 'roc_auc', '-v7.3');
        end
    end
end

end