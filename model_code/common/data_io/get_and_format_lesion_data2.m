function cfg = get_and_format_lesion_data2(beh_table, id_col, beh_col, cfg)

% Load behavioral csv file
if ~istable(beh_table)
    if ~contains(beh_table, '.xls')
        beh_table = readtable(beh_table, 'NumHeaderLines', 0, 'Delimiter', ',');
    else
        beh_table = readtable(beh_table);
    end
end

% Get study IDs and filepaths
sub_ids = string(beh_table.(id_col));

% Get IDs that have behavioral data
beh_data = beh_table.(beh_col);
has_data = ~isnan(beh_data);
beh_data = beh_data(has_data);
sub_ids = sub_ids(has_data);
inc_inds = find(has_data);
has_data = has_data(has_data);

% Verify that lesion images exist for all patients with behavioral data
lesion_check = zeros(length(sub_ids), 1);
for i = 1:length(sub_ids)
    lesion_check(i) = isfile(fullfile(sub_ids(i)));
end
missing_lesion = find(lesion_check == 0);
if ~isempty(missing_lesion)
    warning([num2str(numel(missing_lesion)) ' patients do not have lesion files in the lesion directory. These patients will be excluded.']);
    sub_ids = sub_ids(lesion_check==1);
    has_data(lesion_check==0)=0;
    beh_data = beh_data(lesion_check==1);
    disp(['After excluding patients with missing lesion files, N=' num2str(length(sub_ids)) ' patients will be included'])
else
    disp(['All subjects with behavioral data have lesion files, N=' num2str(length(sub_ids)) ' patients will be included'])
end

% Loop over subjects
cfg.X = zeros(length(sub_ids), length(cfg.mask_inds)); % preallocate matrix
cfg.lvol = zeros(length(sub_ids), 1); % preallocate array

for i = 1:length(sub_ids)
    
    % Load lesion
    if isfile(fullfile(sub_ids(i)))
       my_lesion = niftiread(fullfile(sub_ids(i)));
    end

    % Extract in-mask voxels and fill matrix row
    cfg.X(i,:) = my_lesion(cfg.mask_inds);
    cfg.lvol(i,1) = sum(my_lesion(:));

end

% Define additional outputs
cfg.Y = beh_data;
cfg.include_subs = sub_ids;
cfg.include_inds = inc_inds(has_data);

end