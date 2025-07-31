function [model_results] = run_perm_mass_uni_corr(X, Y, cfg, model_results)
    
% Permutation tests for mass univariate correlations
% Joseph Griffis 2024

% Get permuted Y
if ~isempty(cfg.confounds)
    mdl = fitlm(cfg.confounds, Y);
    Y = mdl.Residuals.Raw;
end

n_perm = cfg.perm.n_perm;
perm_Y = zeros(length(Y), n_perm);        
for i = 1:n_perm
    perm_Y(:,i) = Y(randperm(length(Y))');
end    

% Preallocate permutation correlation matrix
t_perm = zeros(size(X,2), n_perm);

% Run permutations
if ~isfield(cfg, 'cor_type')
    cor_type = 'Pearson';
else
    cor_type = cfg.cor_type;
end
if isempty(cfg.confounds)
    if cfg.parallel == 1
        parfor i = 1:n_perm
            disp(['Running permutation analysis: ' num2str(i)]);
            r = corr(X, perm_Y(:,i), 'type', cor_type);
            t_perm(:,i) = r ./ sqrt((1-r.^2) ./ (length(Y)-2));
        end
    else
        for i = 1:n_perm
            disp(['Running permutation analysis: ' num2str(i)]);
            r = corr(X, perm_Y(:,i), 'type', cor_type);
            t_perm(:,i) = r ./ sqrt((1-r.^2) ./ (length(Y)-2));
        end
    end
else
    confounds = cfg.confounds;
    if cfg.parallel == 1
        parfor i = 1:n_perm
            disp(['Running permutation analysis: ' num2str(i)]);
            r = partialcorr(X, perm_Y(:,i), confounds, 'type', cor_type);
            t_perm(:,i) = r .* sqrt((length(Y)-2-size(confounds,2)) ./ (1-r.^2));
        end
    else
        for i = 1:n_perm
            disp(['Running permutation analysis: ' num2str(i)]);
            r = partialcorr(X, perm_Y(:,i), confounds, 'type', cor_type);
            t_perm(:,i) = r .* sqrt((length(Y)-2-size(confounds,2)) ./ (1-r.^2));
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

