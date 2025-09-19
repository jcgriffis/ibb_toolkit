function cfg = get_and_format_array_data2(beh_table, id_col, beh_col, cfg)

% Load behavioral csv file
if ~istable(beh_table)
    beh_table = readtable(beh_table, 'NumHeaderLines', 0, 'Delimiter', ',');
end

% Get study IDs
sub_ids = string(beh_table.(id_col));

% Get IDs that have behavioral data
beh_data = beh_table.(beh_col);
has_data = ~isnan(beh_data);
beh_data = beh_data(has_data);
sub_ids = sub_ids(has_data);
inc_inds = find(has_data);
has_data = has_data(has_data);

% Verify that lesion images exist for all patients with behavioral data
data_check = zeros(length(sub_ids), 1);
for i = 1:length(sub_ids)
    data_check(i) = isfile(fullfile(sub_ids(i)));
end
missing_data = find(data_check == 0);
if ~isempty(missing_data)
    warning([num2str(numel(missing_data)) ' patients do not have data files in the lesion directory. These patients will be excluded.']);
    sub_ids = sub_ids(data_check==1);
    beh_data = beh_data(data_check==1);
    has_data(data_check==0) = 0;
    disp(['After excluding patients with missing data files, N=' num2str(length(sub_ids)) ' patients will be included'])
else
    disp(['All subjects with behavioral data have data files, N=' num2str(length(sub_ids)) ' patients will be included'])
end

% Loop over subjects
for i = 1:length(sub_ids)
    
    % Load lesion
    if isfile(fullfile(sub_ids(i)))
       my_data = readtable(fullfile(sub_ids(i)));
    end

    % Extract in-mask voxels and fill matrix row
    cfg.X(i,:) = my_data{:,:};

end

% Define outputs
cfg.Y = beh_data;
cfg.include_subs = sub_ids;
cfg.include_inds = inc_inds(has_data);

end