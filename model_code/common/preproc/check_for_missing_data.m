function [has_data, X_trim, Y, cfg] = check_for_missing_data(X_trim, Y, cfg)

% Checks predictor matrix for missing data after trimming, and removes any 
% subjects with no predictor data remaining after trimming
% Joseph Griffis 2023

% Check to see if any subjects have no viable lesion data after trimming X
has_data = nansum(X_trim, 2) > 0;

% Remove any subjects without viable lesion data from X and Y
if ~isempty(has_data(has_data==0))
    disp([num2str(numel(has_data(has_data==0))) ' subjects have no predictor data after trimming...removing from analysis']);
    
    % Remove subjects with no predictor observatiosn from X and Y
    X_trim = X_trim(has_data~=0,:);
    Y = Y(has_data~=0);
    
    % If there is a stratification variable, remove them from it too
    if isfield(cfg, 'strat_var') && ~isempty(cfg.strat_var)
        cfg.strat_var = cfg.strat_var(has_data);
    end
    
    % If confound variables are included, remove subjects from them too
    if isfield(cfg, 'confounds') && ~isempty(cfg.confounds)
        cfg.confounds = cfg.confounds(has_data,:);
    end
    
    % Same for lesion volume
    if isfield(cfg, 'lvol') && ~isempty(cfg.lvol)
        cfg.lvol = cfg.lvol(has_data,:);
    end
    
else
    disp(['All subjects still have predictor data after trimming, proceding with analysis for all subjects (N=' num2str(size(X_trim,1))]);
end

end