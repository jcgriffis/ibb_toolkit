function write_connectome_weights_ulsm(model_results)

% Get thresholds
fwe_thresh = strsplit(num2str(model_results.cfg.fwe_thresh), '.');             
fwe_thresh = fwe_thresh{2};
fdr_thresh = strsplit(num2str(model_results.cfg.fdr_thresh), '.');             
fdr_thresh = fdr_thresh{2};    
unc_thresh = strsplit(num2str(model_results.cfg.unc_thresh), '.');
unc_thresh = unc_thresh{2};

% Unthresholded coeff weight map 
stat = model_results.coeff;
out_name = 'corr_map_unthresh';
write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);

% Unthresholded t-statistic map
stat = model_results.tstat;
out_name = 'tstat_map_unthresh';   
write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);

% FWE results for theoretical p-values
if isfield(model_results, 'coeff_vfwe_pvals')
    stat = model_results.tstat .* (model_results.coeff_vfwe_pvals < model_results.cfg.fwe_thresh);
    out_name = ['tstat_vfwe' fwe_thresh '.nii'];
    write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);

    stat = model_results.coeff .* (model_results.coeff_vfwe_pvals < model_results.cfg.fwe_thresh);
    out_name = ['coeff_vfwe' fwe_thresh '.nii'];
    write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);    
end

% FDR results for theoretical p-values
if isfield(model_results, 'coeff_vfdr_pvals')
    stat = model_results.tstat .* (model_results.coeff_vfdr_pvals < model_results.cfg.fdr_thresh);
    out_name = ['tstat_vfdr' fdr_thresh '.nii'];
    write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);

    stat = model_results.coeff .* (model_results.coeff_vfdr_pvals < model_results.cfg.fdr_thresh);
    out_name = ['coeff_vfdr' fdr_thresh '.nii'];
    write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);    
end

% Permutation results
if isfield(model_results.cfg, 'perm') && model_results.cfg.perm.write_perm_images == 1
    cd(model_results.cfg.out_dir);
    if ~isfolder('Permutation_Result_Images')
        mkdir('Permutation_Result_Images')
    end
    cd(fullfile(model_results.cfg.out_dir, 'Permutation_Result_Images'));
    % Uncorrected voxel p-values
    if model_results.cfg.perm.write_uncorrected_images == 1 && isfield(model_results.perm, 'p_coeff')
        stat = model_results.tstat .* (model_results.perm.p_coeff < model_results.cfg.unc_thresh);
        out_name = ['tstat_uncp' unc_thresh '.nii'];
        write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);

        stat = model_results.coeff .* (model_results.perm.p_coeff < model_results.cfg.unc_thresh);
        out_name = ['coeff_uncp' unc_thresh '.nii'];       
        write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);        
    end

    % FWE voxel p-values
    if isfield(model_results.perm, 'fwep_coeff')
        stat = model_results.tstat .* (model_results.perm.fwep_coeff < model_results.cfg.fwe_thresh);
        out_name = ['tstat_vfwe' fwe_thresh '.nii'];
        write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);
        
        stat = model_results.coeff .* (model_results.perm.fwep_coeff < model_results.cfg.fwe_thresh);
        out_name = ['coeff_vfwe' fwe_thresh '.nii'];
        write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);      
    end
    
    % FWE voxel p-values
    if isfield(model_results.perm, 'fdrp_coeff')
        stat = model_results.tstat .* (model_results.perm.fdrp_coeff < model_results.cfg.fwe_thresh);
        out_name = ['tstat_vfdr' fdr_thresh '.nii'];
        write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);

        stat = model_results.coeff .* (model_results.perm.fdrp_coeff < model_results.cfg.fwe_thresh);
        out_name = ['coeff_vfdr' fdr_thresh '.nii'];
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
                coeff_thresh = model_results.tstat;
                coeff_thresh(abs(coeff_thresh) <= crit) = 0;
                my_coeff = model_results.coeff .* (coeff_thresh ~= 0);
                stat = my_coeff;
                out_name = ['coeff_fwe' fwe_thresh '_v' crit_vals{i}];
                write_edge_file(model_results.cfg.parcel_table, model_results.cfg.tu_mask, model_results.X_ind, stat, model_results.cfg.direction, out_name);               
            end
        end
    end
end
