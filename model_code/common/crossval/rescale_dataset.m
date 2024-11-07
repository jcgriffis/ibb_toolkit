function [cfg, x_train, x_test, y_train, y_test] = rescale_dataset(cfg, x_train, x_test, y_train, y_test)

% Standardization method
if ~isfield(cfg, 'standardize_method')
    cfg.standardize_method = 'zscore';
    cfg.standardize_type = 'std';
end

% Apply standardization if specified
if cfg.standardize == 1  % Convert X to z-scores
    miss_col = sum(x_train,1)==0;
    [x_train, cfg.Cx, cfg.Sx] = normalize(x_train, cfg.standardize_method, cfg.standardize_type);
    x_train(:,miss_col)=0; % Since normalize will cause columns to become NaN if they are all 0
    x_test = normalize(x_test, 'center', cfg.Cx, 'scale', cfg.Sx);
    x_test(:,miss_col)=0; % Since normalize will cause columns to become NaN if they have SD = 0 in train set
elseif cfg.standardize == 2 % Convert Y to z-scores 
    if cfg.cat_Y == 1
        cfg.standardize = 0;
        disp('Unable to standardize categorical Y, setting standardize flag to 0');    
    else
        [y_train, cfg.Cy, cfg.Sy] = normalize(y_train, cfg.standardize_method, cfg.standardize_type);
        y_test = normalize(y_test, 'center', cfg.Cy, 'scale', cfg.Sy);
    end
elseif cfg.standardize == 3 % Convert X and Y to z-scores
    miss_col = sum(x_train,1)==0;
    [x_train, cfg.Cx, cfg.Sx] = normalize(x_train, cfg.standardize_method, cfg.standardize_type);
    x_train(:,miss_col)=0; % Since normalize will cause columns to become NaN if they are all 0
    x_test = normalize(x_test, 'center', cfg.Cx, 'scale', cfg.Sx);
    x_test(:,miss_col)=0; % Since normalize will cause columns to become NaN if they have SD = 0 in train set
    if cfg.cat_Y == 0
        disp('Unable to standardize categorical Y, setting standardize flag to 1');
        cfg.standardize = 1;
    else
        [y_train, cfg.Cy, cfg.Sy] = normalize(y_train, cfg.standardize_method, cfg.standardize_type);
        y_test = normalize(y_test, 'center', cfg.Cy, 'scale', cfg.Sy);      
    end
end