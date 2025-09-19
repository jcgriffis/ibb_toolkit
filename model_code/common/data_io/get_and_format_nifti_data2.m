function cfg = get_and_format_nifti_data2(beh_table, id_col, beh_col, cfg)

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

% Verify that NIFTI images exist for all patients with behavioral data
nii_check = zeros(length(sub_ids), 1);
for i = 1:length(sub_ids)
    nii_check(i) = isfile(fullfile(sub_ids(i)));
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

for i = 1:length(sub_ids)
    
    % Load NIFTI
    if isfile(fullfile(sub_ids(i)))
       my_nii = niftiread(fullfile(sub_ids(i)));
    end

    % Extract in-mask voxels and fill matrix row
    cfg.X(i,:) = my_nii(cfg.mask_inds);

end

% Define additional outputs
cfg.Y = beh_data;
cfg.include_subs = sub_ids;
cfg.include_inds = inc_inds(has_data);

end
    