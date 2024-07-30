function write_mlsm_weight_maps(nifti_out, nifti_hdr, model_results)

% Get thresholds
disp(model_results.cfg.fwe_thresh);
fwe_thresh = strsplit(num2str(model_results.cfg.fwe_thresh), '.');             
fwe_thresh = fwe_thresh{2};
fdr_thresh = strsplit(num2str(model_results.cfg.fdr_thresh), '.');             
fdr_thresh = fdr_thresh{2};    
unc_thresh = strsplit(num2str(model_results.cfg.unc_thresh), '.');
unc_thresh = unc_thresh{2};

% Unthresholded coeff weight map 
nifti_out(:) = 0;
nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.coeff;
niftiwrite(nifti_out, 'coeff_map_unthresh', nifti_hdr);     

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
           

nifti_out(:) = 0;
nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = get_summary(model_results.coeff);
niftiwrite(nifti_out, [prefix '_unthresh'], nifti_hdr);    

% Jack-knife thresholded maps
if isfield(model_results.cfg, 'jack') && model_results.cfg.jack.write_jack_images == 1
    cd(model_results.cfg.out_dir);
    if ~isfolder('Jackknife_Result_Images')
        mkdir('Jackknife_Result_Images')
    end
    cd(fullfile(model_results.cfg.out_dir, 'Jackknife_Result_Images'));
    
    % FWE correction
    if model_results.cfg.jack.get_pvals == 1
       
       % FWE correction
        if model_results.cfg.jack.coeff_fwe == 1
           nifti_out(:) = 0;
           nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = get_summary(model_results.coeff) .* (model_results.jack.fwep_coeff <= model_results.cfg.fwe_thresh);
           niftiwrite(nifti_out, [prefix '_fwep' fwe_thresh], nifti_hdr);      
        end
    
        % FDR correction
        if model_results.cfg.jack.coeff_fdr == 1
           nifti_out(:) = 0;
           nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = get_summary(model_results.coeff) .* (model_results.jack.fdrp_coeff <= model_results.cfg.fdr_thresh);
           niftiwrite(nifti_out, [prefix '_fdrp' fdr_thresh], nifti_hdr);             
        end
        
        % No correction
        if model_results.cfg.jack.write_uncorrected_images == 1
           nifti_out(:) = 0;
           nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = get_summary(model_results.coeff) .* (model_results.jack.p_coeff <= model_results.cfg.unc_thresh);
           niftiwrite(nifti_out, [prefix '_uncp' unc_thresh], nifti_hdr);      
        end
        
        % T-statistic maps

        % Unthresholded t-statistic image
        nifti_out(:) = 0;
        nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.jack.t_coeff;
        niftiwrite(nifti_out, 'tstat_unthresh', nifti_hdr); 

        % FWE correction
        if model_results.cfg.jack.coeff_fwe == 1
           nifti_out(:) = 0;
           nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.jack.t_coeff .* (model_results.jack.fwep_coeff <= model_results.cfg.fwe_thresh);
           niftiwrite(nifti_out, ['tstat_fwep' fwe_thresh], nifti_hdr);      
        end
    
        % FDR correction
        if model_results.cfg.jack.coeff_fdr == 1
           nifti_out(:) = 0;
           nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.jack.t_coeff .* (model_results.jack.fdrp_coeff <= model_results.cfg.fdr_thresh);
           niftiwrite(nifti_out, ['tstat_fdrp' fdr_thresh], nifti_hdr);             
        end
        
        % No correction
        if model_results.cfg.jack.write_uncorrected_images == 1
           nifti_out(:) = 0;
           nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.jack.t_coeff .* (model_results.jack.p_coeff <= model_results.cfg.unc_thresh);
           niftiwrite(nifti_out, ['tstat_uncp' unc_thresh], nifti_hdr);      
        end
    end
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
       
       % FWE correction
        if model_results.cfg.boot.coeff_fwe == 1
           nifti_out(:) = 0;
           nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = get_summary(model_results.coeff) .* (model_results.boot.fwep_coeff <= model_results.cfg.fwe_thresh);
           niftiwrite(nifti_out, [prefix '_fwep' fwe_thresh], nifti_hdr);      
        end
    
        % FDR correction
        if model_results.cfg.boot.coeff_fdr == 1
           nifti_out(:) = 0;
           nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = get_summary(model_results.coeff) .* (model_results.boot.fdrp_coeff <= model_results.cfg.fdr_thresh);
           niftiwrite(nifti_out, [prefix '_fdrp' fdr_thresh], nifti_hdr);             
        end
        
        % No correction
        if model_results.cfg.boot.write_uncorrected_images == 1
           nifti_out(:) = 0;
           nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = get_summary(model_results.coeff) .* (model_results.boot.p_coeff <= model_results.cfg.unc_thresh);
           niftiwrite(nifti_out, [prefix '_uncp' unc_thresh], nifti_hdr);      
        end

       % Z-statistic maps

       % Unthresholded z-statistic image
       nifti_out(:) = 0;
       nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.boot.z_coeff;
       niftiwrite(nifti_out, 'zstat_unthresh', nifti_hdr); 
        
       % FWE correction
        if model_results.cfg.boot.coeff_fwe == 1
           nifti_out(:) = 0;
           nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.boot.z_coeff .* (model_results.boot.fwep_coeff <= model_results.cfg.fwe_thresh);
           niftiwrite(nifti_out, ['zstat_fwep' fwe_thresh], nifti_hdr);      
        end
    
        % FDR correction
        if model_results.cfg.boot.coeff_fdr == 1
           nifti_out(:) = 0;
           nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.boot.z_coeff .* (model_results.boot.fdrp_coeff <= model_results.cfg.fdr_thresh);
           niftiwrite(nifti_out, ['zstat_fdrp' fdr_thresh], nifti_hdr);             
        end
        
        % No correction
        if model_results.cfg.boot.write_uncorrected_images == 1
           nifti_out(:) = 0;
           nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = model_results.boot.z_coeff  .* (model_results.boot.p_coeff <= model_results.cfg.unc_thresh);
           niftiwrite(nifti_out, ['zstat_uncp' unc_thresh], nifti_hdr);      
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
        nifti_out(:) = 0;
        nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = get_summary(model_results.coeff) .* (model_results.perm.p_coeff <= model_results.cfg.unc_thresh);
        niftiwrite(nifti_out, [prefix '_uncp' unc_thresh '.nii'], nifti_hdr);        
    end

    % FWE voxel p-values
    if isfield(model_results.perm, 'fwep_coeff')
        nifti_out(:) = 0;
        nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = get_summary(model_results.coeff) .* (model_results.perm.fwep_coeff <= model_results.cfg.fwe_thresh);
        niftiwrite(nifti_out, [prefix '_vfwe' fwe_thresh '.nii'], nifti_hdr);        
    end
    
    % FWE voxel p-values
    if isfield(model_results.perm, 'fdrp_coeff')
        nifti_out(:) = 0;
        nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = get_summary(model_results.coeff) .* (model_results.perm.fdrp_coeff <= model_results.cfg.fwe_thresh);
        niftiwrite(nifti_out, [prefix '_vfdr' fdr_thresh '.nii'], nifti_hdr);        
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
                nifti_out(:) = 0;
                nifti_out(model_results.cfg.mask_inds(model_results.X_ind)) = my_coeff;
                niftiwrite(nifti_out, [prefix '_' crit_vals{i}], nifti_hdr);                
            end
        end
    end
end
