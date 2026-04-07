function [model_results] = run_perm_mass_uni_corr_interaction(X, Y, cfg, model_results)
    
% Permutation tests for mass univariate correlations
% Joseph Griffis 2024

% Preallocate permutation correlation matrix
n_perm = cfg.perm.n_perm;
t_perm = zeros(size(X,2), n_perm);

% Run permutations
if ~isfield(cfg, 'cor_type')
    cor_type = 'Pearson';
else
    cor_type = cfg.cor_type;
end

if cfg.parallel == 1
    confounds = cfg.confounds;
    int_term = confounds(:,cfg.int_term);
    parfor i = 1:size(X,2)
        if strcmp(cor_type, 'Spearman')
          mdl = fitlm(tiedrank([confounds, X(:,i)]), tiedrank(Y));  
          perm_y = mdl.Residuals.Raw(randperm(length(Y))');    
          r = partialcorr(tiedrank(X(:,i)).*tiedrank(int_term), tiedrank(perm_y), tiedrank([X(:,i), confounds]));
          t_perm(:,i) = r .* sqrt((length(Y)-2-(size(confounds,2)+1)) ./ (1-r.^2));          
        else
          mdl = fitlm([confounds, X(:,i)], Y);  
          perm_y = mdl.Residuals.Raw(randperm(length(Y))');    
          r = partialcorr(X(:,i).*int_term, perm_y, [X(:,i), confounds], 'type', cor_type);
          t_perm(:,i) = r .* sqrt((length(Y)-2-(size(confounds,2)+1)) ./ (1-r.^2));               
        end           
    end    
else
    confounds = cfg.confounds;
    int_term = confounds(:,cfg.int_term);
    for i = 1:size(X,2)
        if strcmp(cor_type, 'Spearman')
          mdl = fitlm(tiedrank([confounds, X(:,i)]), tiedrank(Y));  
          perm_y = mdl.Residuals.Raw(randperm(length(Y))');    
          r = partialcorr(tiedrank(X(:,i)).*tiedrank(int_term), tiedrank(perm_y), tiedrank([X(:,i), confounds]));
          t_perm(:,i) = r .* sqrt((length(Y)-2-(size(confounds,2)+1)) ./ (1-r.^2));          
        else
          mdl = fitlm([confounds, X(:,i)], Y);  
          perm_y = mdl.Residuals.Raw(randperm(length(Y))');    
          r = partialcorr(X(:,i).*int_term, perm_y, [X(:,i), confounds], 'type', cor_type);
          t_perm(:,i) = r .* sqrt((length(Y)-2-(size(confounds,2)+1)) ./ (1-r.^2));               
        end        
    end              
end

% Get permutation results
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