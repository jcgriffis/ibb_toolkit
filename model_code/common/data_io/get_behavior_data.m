function [beh_data, include_ids, include_inds, include_beh] = get_behavior_data(beh_table, beh, subject_list)

% This function is a convenience function for subsetting out behavioral data from spreadsheets 
% and identifying the corresponding subject IDs and indices in the pre-generated lesion matrix
% Inputs:
% 
% beh_table - table containing behavioral data
% beh - name of the column corresponding to the variable of interest
% subject_list - list of subject IDs loaded from dataset file
%
% Outputs:
%
% beh_data - column of scores for subjects without missing data
% include_ids - patient IDs for included subjects in beh_data
% include_inds - row indices for included subjects in pre-generated lesion matrix (or other pre-generated imaging datasets)
% include_beh - row indices for included subject in behavioral data table

% Joseph Griffis 2024

% Get column from table
if ~isstring(beh)
    beh = string(beh);
end

% Get intersection of behavioral data subject list with registry list
beh_id = beh_table{:,1};
sub_id = str2double(subject_list);
[~, ia, ib] = intersect(beh_id, sub_id);

% Subset behavioral data to include only subjects with lesions in the registry
beh_table = beh_table(ia, :);
subject_list = subject_list(ib);

% Find subjects with data
my_col = beh_table.(beh);
has_data = find(~isnan(my_col));

% Get included indices
include_inds = ib(has_data);
include_ids = subject_list(has_data);
include_beh = ia(has_data);

% Subset data
beh_data = my_col(has_data);

end
