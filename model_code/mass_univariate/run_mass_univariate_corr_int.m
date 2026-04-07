function [model_results] = run_mass_univariate_corr_interaction(X, Y, cfg, model_results)

% Run mass-univariate correlation analyses
% Joseph Griffis 2024

% Run mass univariate correlations
if ~isfield(cfg, 'cor_type')
    cor_type = 'Pearson';
else
    cor_type = cfg.cor_type;
end
coeff = zeros(size(X,2),1);
coeff_pvals = zeros(size(X,2),1);

disp('Confounds specified...running partial correlations')
if cfg.parallel == 1
    confounds = cfg.confounds;
    int_term = confounds(:,cfg.int_term);
    parfor i = 1:size(X,2)
        if strcmp(cor_type, 'Spearman')
            [coeff(i,1), coeff_pvals(i,1)] = partialcorr(tiedrank(X(:,i)).*tiedrank(int_term), tiedrank(Y), tiedrank([X(:,i), confounds]));
        else
            [coeff(i,1), coeff_pvals(i,1)] = partialcorr(X(:,i).*int_term, Y, [X(:,i), confounds]);
        end                    
    end    
    model_results.coeff = coeff;
    model_results.coeff_pvals = coeff_pvals;
    clear coeff coeff_pvals
else
    confounds = cfg.confounds;
    int_term = confounds(:,cfg.int_term);            
    for i = 1:size(X,2)
        if strcmp(cor_type, 'Spearman')
            [coeff(i,1), coeff_pvals(i,1)] = partialcorr(tiedrank(X(:,i)).*tiedrank(int_term), tiedrank(Y), tiedrank([X(:,i), confounds]));
        else
            [coeff(i,1), coeff_pvals(i,1)] = partialcorr(X(:,i).*int_term, Y, [X(:,i), confounds]);
        end    
    end
end

% Correct voxel-level p-values if indicated
model_results.coeff_vfwe_pvals = bonf_holm(model_results.coeff_pvals);
model_results.coeff_vfdr_pvals = mafdr(model_results.coeff_pvals, 'BHFDR', true);

% Get t-statistics
if isempty(cfg.confounds)
    model_results.tstat = model_results.coeff ./ sqrt((1-model_results.coeff.^2) ./ (length(Y)-2));
else
    model_results.tstat = model_results.coeff .* sqrt((length(Y)-2-size(cfg.confounds,2)) ./ (1-model_results.coeff.^2));
end

if isempty(cfg.confounds)
    % Generate CPM model
    model_results.obs_y = Y;
    cpm_train = X*model_results.coeff;
    if numel(unique(Y)) > 2
        mdl = fitlm(cpm_train, Y);
    else
        mdl = fitlm(Y, cpm_train, 'CategoricalVars', 'x1');
    end
    model_results.cpm_mdl = mdl;
end

end