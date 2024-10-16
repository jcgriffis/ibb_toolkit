function write_mass_uni_weight_maps(nifti_out, nifti_hdr, model_results)

if ~strcmp(model_results.cfg.model_spec, 'prop_sub')

    % Get thresholds
    fwe_thresh = strsplit(num2str(model_results.cfg.fwe_thresh), '.');             
    fwe_thresh = fwe_thresh{2};
    fdr_thresh = strsplit(num2str(model_results.cfg.fdr_thresh), '.');             
    fdr_thresh = fdr_thresh{2};    
    unc_thresh = strsplit(num2str(model_results.cfg.unc_thresh), '.');
    unc_thresh = unc_thresh{2};

    % Unthresholded coeff weight map 
    if ~strcmp(model_results.cfg.model_spec, 'munilr') && ~ strcmp(model_results.cfg.model_spec, 'ttest') && ~strcmp(model_results.cfg.model_spec, 'muniolsr')
        nifti_out(:) = 0;
        nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.coeff;
        niftiwrite(nifti_out, 'corr_map_unthresh', nifti_hdr);     
    end
    
    % Unthresholded t-statistic map
    nifti_out(:) = 0;
    if strcmp(model_results.cfg.model_spec, 'bmunz')
        model_results.tstat = model_results.coeff;
        stat_name = 'bmstat';
    elseif strcmp(model_results.cfg.model_spec, 'ttest')
        stat_name = 'tstat';
    else
        stat_name = 'tstat';
    end
    nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.tstat;
    niftiwrite(nifti_out, [stat_name '_map_unthresh'], nifti_hdr);    
    
    % FWE results for theoretical p-values
    if isfield(model_results, 'coeff_vfwe_pvals')
        nifti_out(:) = 0;
        nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.tstat .* (model_results.coeff_vfwe_pvals < model_results.cfg.fwe_thresh);
        niftiwrite(nifti_out, [stat_name '_vfwe' fwe_thresh '.nii'], nifti_hdr);
        if ~strcmp(model_results.cfg.model_spec, 'munilr') && ~strcmp(model_results.cfg.model_spec, 'muniolsr') && ~strcmp(model_results.cfg.model_spec, 'bmunz') && ~strcmp(model_results.cfg.model_spec, 'ttest') 
            nifti_out(:) = 0;
            nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.coeff .* (model_results.coeff_vfwe_pvals < model_results.cfg.fwe_thresh);
            niftiwrite(nifti_out, ['coeff_vfwe' fwe_thresh '.nii'], nifti_hdr);      
        end
    end
    
    % FDR results for theoretical p-values
    if isfield(model_results, 'coeff_vfdr_pvals')
        nifti_out(:) = 0;
        nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.tstat .* (model_results.coeff_vfdr_pvals < model_results.cfg.fdr_thresh);
        niftiwrite(nifti_out, [stat_name '_vfdr' fdr_thresh '.nii'], nifti_hdr);
        if ~strcmp(model_results.cfg.model_spec, 'munilr') && ~strcmp(model_results.cfg.model_spec, 'muniolsr') && ~strcmp(model_results.cfg.model_spec, 'bmunz') && ~strcmp(model_results.cfg.model_spec, 'ttest') 
            nifti_out(:) = 0;
            nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.coeff .* (model_results.coeff_vfdr_pvals < model_results.cfg.fdr_thresh);
            niftiwrite(nifti_out, ['coeff_vfdr' fdr_thresh '.nii'], nifti_hdr);   
        end
    end
    
    % Permutation results
    if isfield(model_results, 'perm')
        if model_results.cfg.perm.write_perm_images==1
            cd(model_results.cfg.out_dir);
            if ~isfolder('Permutation_Result_Images')
                mkdir('Permutation_Result_Images')
            end
            cd(fullfile(model_results.cfg.out_dir, 'Permutation_Result_Images'));
            % Uncorrected voxel p-values
            if model_results.cfg.perm.write_uncorrected_images == 1 && isfield(model_results.perm, 'p_coeff')
                nifti_out(:) = 0;
                nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.tstat .* (model_results.perm.p_coeff < model_results.cfg.unc_thresh);
                niftiwrite(nifti_out, [stat_name '_uncp' unc_thresh '.nii'], nifti_hdr);
                nifti_out(:) = 0;
                if ~strcmp(model_results.cfg.model_spec, 'munilr') && ~strcmp(model_results.cfg.model_spec, 'bmunz') && ~strcmp(model_results.cfg.model_spec, 'ttest')
                    nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.coeff .* (model_results.perm.p_coeff < model_results.cfg.unc_thresh);
                    niftiwrite(nifti_out, ['coeff_uncp' unc_thresh '.nii'], nifti_hdr);        
                end
            end
        
            % FWE voxel p-values
            if isfield(model_results.perm, 'fwep_coeff')
                nifti_out(:) = 0;
                nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.tstat .* (model_results.perm.fwep_coeff < model_results.cfg.fwe_thresh);
                niftiwrite(nifti_out, [stat_name '_vfwe' fwe_thresh '.nii'], nifti_hdr);
                nifti_out(:) = 0;
                if ~strcmp(model_results.cfg.model_spec, 'munilr') && ~strcmp(model_results.cfg.model_spec, 'bmunz') && ~strcmp(model_results.cfg.model_spec, 'ttest')
                    nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.coeff .* (model_results.perm.fwep_coeff < model_results.cfg.fwe_thresh);
                    niftiwrite(nifti_out, ['coeff_vfwe' fwe_thresh '.nii'], nifti_hdr);        
                end
            end
            
            % FWE voxel p-values
            if isfield(model_results.perm, 'fdrp_coeff')
                nifti_out(:) = 0;
                nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.tstat .* (model_results.perm.fdrp_coeff < model_results.cfg.fwe_thresh);
                niftiwrite(nifti_out, [stat_name '_vfdr' fdr_thresh '.nii'], nifti_hdr);
                nifti_out(:) = 0;
                if ~strcmp(model_results.cfg.model_spec, 'munilr') && ~strcmp(model_results.cfg.model_spec, 'muniolsr') && ~strcmp(model_results.cfg.model_spec, 'bmunz') && ~strcmp(model_results.cfg.model_spec, 'ttest') 
                    nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.coeff .* (model_results.perm.fdrp_coeff < model_results.cfg.fwe_thresh);
                    niftiwrite(nifti_out, ['coeff_vfdr' fdr_thresh '.nii'], nifti_hdr);        
                end
            end
        
            % Permutation thresholded images
            if model_results.cfg.perm.write_perm_images == 1
                cd(model_results.cfg.out_dir);
                if ~isfolder('Permutation_Result_Images')
                    mkdir('Permutation_Result_Images')
                end
                if strcmp(model_results.cfg.model_spec, 'munilr') || strcmp(model_results.cfg.model_spec, 'ttest') || strcmp(model_results.cfg.model_spec, 'muniolsr')
                    model_results.coeff = model_results.tstat;
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
                        nifti_out(:) = 0;
                        nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = my_coeff;
                        if strcmp(model_results.cfg.model_spec, 'ttest') || strcmp(model_results.cfg.model_spec, 'munilr') || strcmp(model_results.cfg.model_spec, 'muniolsr')
                            niftiwrite(nifti_out, ['tstat_fwe' fwe_thresh '_v' crit_vals{i}], nifti_hdr);   
                        elseif strcmp(model_results.cfg.model_spec, 'bmunz')
                            niftiwrite(nifti_out, ['bmstat_fwe' fwe_thresh '_v' crit_vals{i}], nifti_hdr);   
                        else
                            niftiwrite(nifti_out, ['coeff_fwe' fwe_thresh '_v' crit_vals{i}], nifti_hdr);
                        end
                    end
                end
            end
        end
    end

else

    % Propotional overlap for group 0
    nifti_out(:) = 0;
    nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.overlap_g0;
    niftiwrite(nifti_out, 'prop_overlap_g0', nifti_hdr);     

    % Propotional overlap for group 1
    nifti_out(:) = 0;
    nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.overlap_g1;
    niftiwrite(nifti_out, 'prop_overlap_g1', nifti_hdr);     

    % Proportional subtraction
    nifti_out(:) = 0;
    nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.g1_minus_g0;
    niftiwrite(nifti_out, 'prop_sub_g1_minus_g0', nifti_hdr);      

    if model_results.cfg.permutation == 1

        % Get thresholds
        fwe_thresh = strsplit(num2str(model_results.cfg.fwe_thresh), '.');             
        fwe_thresh = fwe_thresh{2};

        % Permutation thresholded images
        if model_results.cfg.perm.write_perm_images == 1
            cd(model_results.cfg.out_dir);
            if ~isfolder('Permutation_Result_Images')
                mkdir('Permutation_Result_Images')
            end
            model_results.coeff = model_results.g1_minus_g0;
            cd(fullfile(model_results.cfg.out_dir, 'Permutation_Result_Images'));
            if model_results.cfg.perm.coeff_cfwe == 1
                % Critical values for different v thresholds
                crit_vals = fieldnames(model_results.perm.cfwe);    
                for i = 1:length(crit_vals)
                    crit = model_results.perm.cfwe.(crit_vals{i});                          
                    coeff_thresh = model_results.g1_minus_g0;
                    coeff_thresh(abs(coeff_thresh) <= crit) = 0;
                    my_coeff = model_results.coeff .* (coeff_thresh ~= 0);
                    nifti_out(:) = 0;
                    nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = my_coeff;
                    niftiwrite(nifti_out, ['diff_fwe' fwe_thresh '_v' crit_vals{i}], nifti_hdr);                    
                end
            end
        end
    end
end
