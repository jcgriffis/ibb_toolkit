function cfg = get_and_format_matrix_data(beh_table, id_col, beh_col, registry_flag, data_dir, cfg)

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
inc_inds = find(has_data);
has_data = has_data(has_data);

% Verify that lesion images exist for all patients with behavioral data
matrix_check = zeros(length(sub_ids), 1);
for i = 1:length(sub_ids)
    txtfile_check = isfile(fullfile(data_dir, strcat(sub_ids(i), '.txt')));
    csfvile_check = isfile(fullfile(data_dir, strcat(sub_ids(i), '.csv')));
    if csfvile_check == 0 && txtfile_check == 0
        matrix_check(i) = 0;
    else
        matrix_check(i) = 1;
    end
end
missing_matrix = find(matrix_check == 0);
if ~isempty(missing_matrix)
    warning([num2str(numel(missing_matrix)) ' patients do not have matrix files in the lesion directory. These patients will be excluded.']);
    sub_ids = sub_ids(matrix_check==1);
    beh_data = beh_data(matrix_check==1);
    has_data(matrix_check == 0) = 0;
    disp(['After excluding patients with missing matrix files, N=' num2str(length(sub_ids)) ' patients will be included'])
else
    disp(['All subjects with behavioral data have matrix files, N=' num2str(length(sub_ids)) ' patients will be included'])
end

% Loop over subjects
for i = 1:length(sub_ids)
    
    % Load lesion
    if isfile(fullfile(data_dir, strcat(sub_ids(i), '.txt')))
       my_data = readmatrix(fullfile(data_dir, strcat(sub_ids(i), '.txt')));
    elseif isfile(fullfile(data_dir, strcat(sub_ids(i), '.csv')))
       my_data = readmatrix(fullfile(data_dir, strcat(sub_ids(i), '.csv')));
    end
    
    % Get mask for upper triangular portion of matrix
    tu_mask = triu(ones(size(my_data)), 1);
    cfg.tu_mask = tu_mask;
    cfg.mask_inds = find(tu_mask);

    % Extract in-mask voxels and fill matrix row
    cfg.X(i,:) = my_data(cfg.mask_inds);

end

% Define outputs
cfg.Y = beh_data;
cfg.include_subs = sub_ids;
cfg.include_inds = inc_inds(has_data);

end
    