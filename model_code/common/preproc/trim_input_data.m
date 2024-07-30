function [X_trim, X_inc_ind, Y, cfg] = trim_input_data(cfg)
    
% Remove columns from predictor matrix with insufficient observations
% Joseph Griffis 2023

if ~strcmp(cfg.model_spec, 'olsr')
    % Set any NaN or Inf data to 0
    xnan_check = sum(isnan(cfg.X));
    xinf_check = sum(isinf(cfg.X));
    ynan_check = sum(isnan(cfg.Y));
    yinf_check = sum(isinf(cfg.Y));
    if xnan_check > 0
       warning('NaN values present in predictor matrix, setting NaN values to zero and proceeding with analysis. Check data processing if not expected.');
       cfg.X(isnan(cfg.X))=0;
    end
    if xinf_check > 0
       warning('Inf values present in predictor matrix, setting Inf values to zero and proceeding with analysis. Check data processing if not expected.');
       cfg.X(isinf(cfg.X))=0;
    end
    if ynan_check > 0
       warning('NaN values present in response data, setting NaN values to zero and proceeding with analysis. Check data processing if not expected.');
       cfg.Y(isnan(cfg.Y))=0;
    end
    if yinf_check > 0
       warning('Inf values present in response data, setting Inf values to zero and proceeding with analysis. Check data processing if not expected.');
       cfg.Y(isinf(cfg.Y))=0;
    end

    % Trim predictor data based on min_obs and freq_thresh
    if cfg.trim_X == 1
        [X_trim, X_inc_ind] = trim_X_data(single(cfg.X), cfg.min_obs, cfg.freq_thresh); % get trimmed X_trim matrix and indices for retained cells
    else
        X_trim = single(cfg.X);
        X_inc_ind = 1:length(X_trim);
    end
    Y = cfg.Y; % final response variable
    
    % Remove raw input data from cfg to reduce memory load
    cfg = rmfield(cfg, 'X');
    cfg = rmfield(cfg, 'Y');
else
    % Don't change X since this could break the link between the predictor matrix and the pre-specified model terms
    X_trim = cfg.X;
    Y = cfg.Y;
    X_inc_ind = 1:1:(size(X_trim,2));

    % Remove raw input data from cfg to reduce memory load
    cfg = rmfield(cfg, 'X');
    cfg = rmfield(cfg, 'Y');
end
    
end