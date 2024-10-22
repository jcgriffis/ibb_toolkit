function Y = check_Y_data(X, Y, cfg)

% Check to make sure that Y data are appropriate given the selected modeling strategy
% Joseph Griffis 2024

if cfg.cat_Y == 1
    if isnumeric(Y)
        if min(Y) ~= -1
            disp('Recoding group with minimum value to -1');
            Y(Y == min(Y)) = -1;
        end
        if max(Y) ~= 1
            disp('Recoding group with maximum value to 1');
            Y(Y == max(Y)) = 1;
        end
    elseif iscategorical(Y) || isstring(Y) || ischar(Y)
        warning('Categorical outcome not coded as -1 and 1; trying to recode. Check data to ensure correctness.')
        unique_vals = unique(Y);
        if numel(unique_vals)==2
            temp_Y = zeros(length(Y),1);
            temp_Y(Y==unique_vals(1))=-1;
            temp_Y(Y==unique_vals(2))=1;
            Y = temp_Y; 
            clear temp_Y;
        else
            error('Categorical outcome contains more than 2 categories. Please check the data and try again.');
        end
    else
        error('Categorical outcome specified, but data are not properly formatted. Please code groups as -1 and 1 and try again.');
    end
elseif cfg.cat_Y == 0
    if isnumeric(Y)
        unique_vals_Y = unique(Y);
        unique_vals_X = unique(X(:));
        if length(unique_vals_Y) <= 2
            if strcmp(cfg.model_spec, 'municorr') && length(unique_vals_X) == 2
                error('Modeling approach assumes continuous outcome variable, but outcome variable only has two values. Check data and use classification approach if necessary.');
            elseif ~strcmp(cfg.model_spec, 'municorr')
                error('Modeling approach assumes continuous outcome variable, but outcome variable only has two values. Check data and use classification approach if necessary.');
            end                
        end
        clear unique_vals_Y unique_vals_X
    else
        error('Continuous modeling approach selected, but outcome data are not numeric. Please check data and try again');
    end
end
