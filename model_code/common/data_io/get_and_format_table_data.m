function cfg = get_and_format_table_data(beh_table, id_col, beh_col, pred_cols, cfg)

% Load behavioral csv file
if ~istable(beh_table)
    beh_table = readtable(beh_table);
end

% Get study IDs
sub_ids = beh_table.(id_col);

if registry_flag == 1
    % Force IDs to character
    my_sub_ids = char(num2str(sub_ids));
    
    % Append zeros to IDs less than 1000 
    for i = 1:length(sub_ids)
        if sub_ids(i) < 100
            my_sub_ids((i), 1:2) = '00';
        elseif sub_ids(i) > 100 && sub_ids(i) < 1000
            my_sub_ids((i), 1) = '0';
        end
    end
else 
    my_sub_ids = sub_ids;
end
sub_ids = string(my_sub_ids); % convert from char to string

% Get IDs that have behavioral data
beh_data = beh_table.(beh_col);
has_data = ~isnan(beh_data);
beh_data = beh_data(has_data);
sub_ids = sub_ids(has_data);

% Get predictor data
for i = 1:length(pred_cols)
    pred_data(:,i) = beh_table.(pred_cols);
end

% Define outputs
cfg.X = pred_data(has_data,:);
cfg.X_names = pred_cols;
cfg.Y = beh_data;
cfg.include_subs = sub_ids;

end
    