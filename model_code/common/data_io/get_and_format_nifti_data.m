function cfg = get_and_format_nifti_data(beh_table, id_col, beh_col, registry_flag, nii_dir, cfg)

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
has_data = has_data(has_data);

% Verify that NIFTI images exist for all patients with behavioral data
nii_check = zeros(length(sub_ids), 1);
for i = 1:length(sub_ids)
    nii_gz_check = isfile(fullfile(nii_dir, strcat(sub_ids(i), '.nii.gz')));
    nii_uz_check = isfile(fullfile(nii_dir, strcat(sub_ids(i), '.nii')));
    if nii_uz_check == 0 && nii_gz_check == 0
        nii_check(i) = 0;
    else
        nii_check(i) = 1;
    end
end
missing_nii = find(nii_check == 0);
if ~isempty(missing_nii)
    warning([num2str(numel(missing_nii)) ' patients do not have NIFTI files in the NIFTI directory. These patients will be excluded.']);
    sub_ids = sub_ids(nii_check==1);
    beh_data = beh_data(nii_check==1);
    has_data(nii_check==0) = 0;
    disp(['After excluding patients with missing NIFTI files, N=' num2str(length(sub_ids)) ' patients will be included'])
else
    disp(['All subjects with behavioral data have NIFTI files, N=' num2str(length(sub_ids)) ' patients will be included'])
end

% Loop over subjects
cfg.X = zeros(length(sub_ids), length(cfg.mask_inds)); % preallocate matrix
cfg.lvol = zeros(length(sub_ids), 1); % preallocate array

for i = 1:length(sub_ids)
    
    % Load NIFTI
    if isfile(fullfile(nii_dir, strcat(sub_ids(i), '.nii.gz')))
       my_nii = niftiread(fullfile(nii_dir, strcat(sub_ids(i), '.nii.gz')));
    elseif isfile(fullfile(nii_dir, strcat(sub_ids(i), '.nii')))
       my_nii = niftiread(fullfile(nii_dir, strcat(sub_ids(i), '.nii')));
    end

    % Extract in-mask voxels and fill matrix row
    cfg.X(i,:) = my_nii(cfg.mask_inds);

end

% Define additional outputs
cfg.Y = beh_data;
cfg.include_subs = sub_ids;
cfg.include_inds = find(has_data);

end
    
