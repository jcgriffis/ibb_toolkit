function write_connectome_weights_mlsm(model_results)

% Get thresholds
fwe_thresh = strsplit(num2str(model_results.cfg.fwe_thresh), '.');             
fwe_thresh = fwe_thresh{2};
fdr_thresh = strsplit(num2str(model_results.cfg.fdr_thresh), '.');             
fdr_thresh = fdr_thresh{2};    
unc_thresh = strsplit(num2str(model_results.cfg.unc_thresh), '.');
unc_thresh = unc_thresh{2};

% z-score coefficients
if strcmp(model_results.cfg.map_summary, 'zscore')
    get_summary = @zscore;
    prefix = 'coeff_z_rescaled';
elseif strcmp(model_results.cfg.map_summary, 'max')
    if strcmp(model_results.cfg.direction, 'pos')
        get_summary = @(map) map ./ max(map);
    elseif strcmp(model_results.cfg.direction, 'neg')
        get_summary = @(map) map ./ min(map);
    else
        get_summary = @(map) map ./ max(abs(map));
    end
    prefix = 'coeff_max_rescaled';
end
           
% Check for parcellation table
if ~isfield(model_results.cfg, 'parcel_table') 
    disp('Field parcel_table not found in cfg structure; .node files will not be generated')
    model_results.cfg.parcel_table = [];
elseif isfield(model_results.cfg, 'parcel_table') && isempty(model_results.cfg.parcel_table)
    disp('Field parcel_table not found in cfg structure; .node files will not be generated')
end

% Unthresholded coeff weight map  
write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, model_results.coeff, model_results.cfg.direction, 'coeff_unthresh');

% z-score coefficients
write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, get_summary(model_results.coeff), model_results.cfg.direction, [prefix '_coeff_unthresh']);

% FWE results for theoretical p-values
if isfield(model_results, 'coeff_vfwe_pvals')
    stat = get_summary(model_results.coeff) .* (model_results.coeff_vfwe_pvals <= model_results.cfg.fwe_thresh);
    out_name = [prefix '_coeff_vfwe' fwe_thresh];
    write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);
end

% FDR results for theoretical p-values
if isfield(model_results, 'coeff_vfdr_pvals')
    stat = get_summary(model_results.coeff) .* (model_results.coeff_vfdr_pvals <= model_results.cfg.fdr_thresh);
    out_name = [prefix '_coeff_vfdr' fdr_thresh '.nii'];
    write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);    
end

% Bootstrap thresholded maps
if isfield(model_results.cfg, 'boot') && model_results.cfg.boot.write_boot_images == 1
    cd(model_results.cfg.out_dir);
    if ~isfolder('Bootstrap_Result_Images')
        mkdir('Bootstrap_Result_Images')
    end
    cd(fullfile(model_results.cfg.out_dir, 'Bootstrap_Result_Images'));
    
    % FWE correction
    if model_results.cfg.boot.get_pvals == 1
       
       % Unthresholded z-statistic image
       stat = model_results.boot.z_coeff;
       out_name = 'zstat_unthresh';
       write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);

       % FWE correction
        if model_results.cfg.boot.coeff_fwe == 1
           stat = get_summary(model_results.coeff) .* (model_results.boot.fwep_coeff <= model_results.cfg.fwe_thresh);
           out_name = [prefix '_fwep' fwe_thresh];
           write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);
         
           stat = model_results.boot.z_coeff .* (model_results.boot.fwep_coeff <= model_results.cfg.fwe_thresh);
           out_name = ['zstat_fwep' fwe_thresh];
           write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);        
        end
    
        % FDR correction
        if model_results.cfg.boot.coeff_fdr == 1
           stat = get_summary(model_results.coeff) .* (model_results.boot.fdrp_coeff <= model_results.cfg.fdr_thresh);
           out_name = [prefix '_fdrp' fdr_thresh];
           write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);           
           
           stat = model_results.boot.z_coeff.* (model_results.boot.fdrp_coeff <= model_results.cfg.fdr_thresh);
           out_name = ['zstat_fdrp' fdr_thresh];
           write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);           
        end
        
        % No correction
        if model_results.cfg.boot.write_uncorrected_images == 1
           stat = get_summary(model_results.coeff) .* (model_results.boot.p_coeff <= model_results.cfg.unc_thresh);
           out_name = [prefix '_uncp' unc_thresh];
           write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);           

           stat = model_results.boot.z_coeff .* (model_results.boot.p_coeff <= model_results.cfg.unc_thresh);
           out_name = ['zstat_uncp' unc_thresh];
           write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);                  
        end
    end
end

% Permutation results
if isfield(model_results.cfg, 'perm') && model_results.cfg.perm.write_perm_images == 1
    cd(model_results.cfg.out_dir);
    if ~isfolder('Permutation_Result_Images')
        mkdir('Permutation_Result_Images')
    end

    % Uncorrected voxel p-values
    if model_results.cfg.perm.write_uncorrected_images == 1 && isfield(model_results.perm, 'p_coeff')
        stat = get_summary(model_results.coeff) .* (model_results.perm.p_coeff <= model_results.cfg.unc_thresh);
        out_name = [prefix '_coeff_uncp' unc_thresh '.nii'];
        write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);    
    end

    % FWE voxel p-values
    if isfield(model_results.perm, 'fwep_coeff')
        stat = get_summary(model_results.coeff) .* (model_results.perm.fwep_coeff <= model_results.cfg.fwe_thresh);
        out_name = [prefix '_coeff_vfwe' fwe_thresh '.nii'];     
        write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);        
    end
    
    % FWE voxel p-values
    if isfield(model_results.perm, 'fdrp_coeff')
        stat = get_summary(model_results.coeff) .* (model_results.perm.fdrp_coeff <= model_results.cfg.fwe_thresh);
        out_name = [prefix '_coeff_vfdr' fdr_thresh '.nii'];
        write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);     
    end

    % Permutation thresholded images
    if model_results.cfg.perm.write_perm_images == 1
        cd(model_results.cfg.out_dir);
        if ~isfolder('Permutation_Result_Images')
            mkdir('Permutation_Result_Images')
        end
        cd(fullfile(model_results.cfg.out_dir, 'Permutation_Result_Images'));
        if model_results.cfg.perm.coeff_cfwe == 1
            % Critical values for different v thresholds
            crit_vals = fieldnames(model_results.perm.cfwe);    
            for i = 1:length(crit_vals)
                crit = model_results.perm.cfwe.(crit_vals{i});                          
                coeff_thresh = model_results.coeff;
                coeff_thresh(abs(coeff_thresh) <= crit) = 0;
                my_coeff = get_summary(model_results.coeff) .* (coeff_thresh ~= 0);
                stat = my_coeff;
                out_name = [prefix '_' crit_vals{i}];
                write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);                
            end
        end
    end
end