function cfg = get_and_format_lesion_data(beh_table, id_col, beh_col, registry_flag, lesion_dir, cfg)

% Load behavioral csv file
if ~istable(beh_table)
    if ~contains(beh_table, '.xls')
        beh_table = readtable(beh_table, 'NumHeaderLines', 0, 'Delimiter', ',');
    else
        beh_table = readtable(beh_table);
    end
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
lesion_check = zeros(length(sub_ids), 1);
for i = 1:length(sub_ids)
    nii_gz_check = isfile(fullfile(lesion_dir, strcat(sub_ids(i), '.nii.gz')));
    nii_check = isfile(fullfile(lesion_dir, strcat(sub_ids(i), '.nii')));
    if nii_check == 0 && nii_gz_check == 0
        lesion_check(i) = 0;
    else
        lesion_check(i) = 1;
    end
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
    if isfile(fullfile(lesion_dir, strcat(sub_ids(i), '.nii.gz')))
       my_lesion = niftiread(fullfile(lesion_dir, strcat(sub_ids(i), '.nii.gz')));
    elseif isfile(fullfile(lesion_dir, strcat(sub_ids(i), '.nii')))
       my_lesion = niftiread(fullfile(lesion_dir, strcat(sub_ids(i), '.nii')));
    end

    % Extract in-mask voxels and fill matrix row
    cfg.X(i,:) = my_lesion(cfg.mask_inds);
    cfg.lvol(i,1) = sum(my_lesion(cfg.mask_inds));

end

% Define additional outputs
cfg.Y = beh_data;
cfg.include_subs = sub_ids;
cfg.include_inds = inc_inds(has_data);

end
    
