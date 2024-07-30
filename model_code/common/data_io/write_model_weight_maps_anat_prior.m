function write_model_weight_maps_anat_prior(cfg, model_results, predictor)

% Get brain mask and brain indices
load(cfg.data_file, 'parcel_img');

% Create temp image
temp_img = parcel_img;
temp_img.img(:) = 0;

% Subset out parcels included in analysis
parcel_inds = find(model_results.X_ind <= cfg.parcel_count);

% Assign weights to parcels (unthresholded z-scored Betas)
parcel_betas = zscore(model_results.coeff);
parcel_betas = parcel_betas(parcel_inds);
for i = 1:length(parcel_inds)
    temp_img.img(parcel_img.img == parcel_inds(i)) = parcel_betas(i); 
end
cd(cfg.out_dir);
save_nii(temp_img, 'betas_z_rescaled_unthresh.nii.gz');
temp_img.img(:) = 0;

if cfg.boot_coeff_fwe == 1 && cfg.bootstrap == 1 && cfg.bootstrap_pvals == 1    

    if ~isfolder(fullfile(cfg.out_dir, 'Bootstrap_Result_Images'))
        mkdir(fullfile(cfg.out_dir, 'Bootstrap_Result_Images'));
    end
    cd(fullfile(cfg.out_dir, 'Bootstrap_Result_Images'));

    % Assign weights to parcels (unthresholded z-stats)
    parcel_betas = model_results.boot.z_betas(parcel_inds);
    for i = 1:length(parcel_inds)
        temp_img.img(parcel_img.img == parcel_inds(i)) = parcel_betas(i); 
    end
    save_nii(temp_img, 'betas_zstat_unthresh.nii.gz');
    temp_img.img(:) = 0;

    % Assign weights to parcels (FWEp z-stats)
    parcel_betas = model_results.boot.z_betas .* (model_results.boot.fwep_betas < 0.05);
    for i = 1:length(parcel_inds)
        temp_img.img(parcel_img.img == parcel_inds(i)) = parcel_betas(i); 
    end
    save_nii(temp_img, ['parcel_betas_zstat_fwep' num2str(cfg.fwe_thresh) '.nii.gz']);
    temp_img.img(:) = 0;
end

if cfg.boot_coeff_fdr == 1 && cfg.bootstrap == 1 && cfg.bootstrap_pvals == 1    
    
    if ~isfolder(fullfile(cfg.out_dir, 'Bootstrap_Result_Images'))
        mkdir(fullfile(cfg.out_dir, 'Bootstrap_Result_Images'));
    end
    cd(fullfile(cfg.out_dir, 'Bootstrap_Result_Images'));

    % Assign weights to parcels (FDRp z-stats)
    parcel_betas = model_results.boot.z_betas .* (model_results.boot.fdrp_betas < cfg.fdr_thresh);
    for i = 1:length(parcel_inds)
        temp_img.img(parcel_img.img == parcel_inds(i)) = parcel_betas(i); 
    end
    save_nii(temp_img, ['parcel_betas_zstat_fdrp' num2str(cfg.fdr_thresh) '.nii.gz']);
    temp_img.img(:) = 0;
end

if cfg.permutation == 1 && cfg.perm_coeff_cfwe == 1

    if ~isfolder(fullfile(cfg.out_dir, 'Permutation_Result_Images'))
        mkdir(fullfile(cfg.out_dir, 'Permutation_Result_Images'));
    end
    cd(fullfile(cfg.out_dir, 'Permutation_Result_Images'));   

    crit_vals = fieldnames(model_results.perm.beta_cfwe);  
    for i = 1:length(crit_vals)
        
        temp_img.img(:) = 0;

        % Get critical value
        crit = model_results.perm.beta_cfwe.(crit_vals{i});    
    
        % Threshold weights
        parcel_betas = model_results.coeff ./ max(abs(model_results.coeff));
        if strcmp(model_results.cfg.direction, 'pos')
            parcel_betas(parcel_betas < crit) = 0;
        elseif strcmp(model_results.cfg.direction, 'neg')
             parcel_betas(parcel_betas > crit) = 0;
        else
             parcel_betas(abs(parcel_betas) < crit) = 0;
        end
        
        % Assign weights to parcels
        for j = 1:length(parcel_inds)
            temp_img.img(parcel_img.img == parcel_inds(j)) = parcel_betas(j); 
        end
        save_nii(temp_img, ['parcel_betas_max_rescaled_cfwe' num2str(cfg.fwe_thresh) '_' crit_vals{i} '.nii.gz']);

    end
end

if strcmp(predictor, 'anat_prior')

    % Tract figures
    tract_inds = find(model_results.X_ind >= cfg.parcel_count);
    tract_betas = zscore(model_results.coeff);
    tract_betas = tract_betas(tract_inds);
    tract_name = cfg.predictor_names(model_results.X_ind(tract_inds));
    
    % Make figures
    scrsz = get(groot,'ScreenSize');
    
    % Sorted betas (Z-rescaled)
    fig = figure('Position',[1 scrsz(4)/3 scrsz(3)/4 scrsz(4)/2.5]);
    [sorted_betas, sort_ind] = sort(tract_betas, 'descend');
    bar(sorted_betas);
    set(gca, 'Xtick', 1:1:length(tract_betas))
    set(gca, 'XtickLabel', tract_name(sort_ind));
    set(gca, 'XtickLabelRotation', 90);
    set(gca, 'TickLabelInterpreter', 'none');
    title('Tract Weights')
    ylabel('Z-Score');
    set(gcf, 'color', 'w');
    grid on;
    saveas(fig, 'sorted_z_rescaled_tract_betas.fig');
    close all;
    
    % Sorted betas (FWEp)
    fig = figure('Position',[1 scrsz(4)/3 scrsz(3)/4 scrsz(4)/2.5]);
    [sorted_betas, sort_ind] = sort(tract_betas, 'descend');
    fwep_vals = model_results.boot.fwep_betas(tract_inds);
    fwep_vals = fwep_vals(sort_ind);
    sorted_betas = sorted_betas(fwep_vals < cfg.fwe_thresh);
    bar(sorted_betas);
    set(gca, 'Xtick', 1:1:length(sorted_betas))
    set(gca, 'XtickLabel', tract_name(sort_ind(fwep_vals < cfg.fwe_thresh)));
    set(gca, 'XtickLabelRotation', 90);
    set(gca, 'TickLabelInterpreter', 'none');
    title(['Tract Weights (FWEp < ' num2str(cfg.fwe_thresh) ')'])
    ylabel('Z-Scored Weights');
    set(gcf, 'color', 'w');
    grid on;
    saveas(fig, 'sorted_z_rescaled_tract_betas_fwep.fig');
    close all;
    
    % Sorted betas (FDRp)
    fig = figure('Position',[1 scrsz(4)/3 scrsz(3)/4 scrsz(4)/2.5]);
    [sorted_betas, sort_ind] = sort(tract_betas, 'descend');
    fdrp_vals = model_results.boot.fdrp_betas(tract_inds);
    fdrp_vals = fdrp_vals(sort_ind);
    sorted_betas = sorted_betas(fdrp_vals < cfg.fdr_thresh);
    bar(sorted_betas);
    set(gca, 'Xtick', 1:1:length(sorted_betas))
    set(gca, 'XtickLabel', tract_name(sort_ind(fdrp_vals < cfg.fdr_thresh)));
    set(gca, 'XtickLabelRotation', 90);
    set(gca, 'TickLabelInterpreter', 'none');
    title(['Tract Weights (FDRp < ' num2str(cfg.fdr_thresh) ')'])
    ylabel('Z-Scored Weights');
    set(gcf, 'color', 'w');
    grid on;
    saveas(fig, 'sorted_z_rescaled_tract_betas_fdrp.fig');
    close all;
    
    % Comparison of Parcel vs. Tract Weights
    fig = figure('Position',[1 scrsz(4)/3 scrsz(3)/13 scrsz(4)/4]);
    boxchart([ones(length(parcel_inds),1);ones(length(tract_inds),1).*2], zscore(model_results.coeff)')
    grid on;
    set(gcf, 'color', 'w');
    ylabel('Z-Scored Weights');
    set(gca, 'Xtick', 1:1:2);
    set(gca, 'XtickLabel', {'Parcels', 'Tracts'})
    title('Parcel vs. Tract Weight Distributions', 'FontSize', 8);
    saveas(fig, 'parcel_vs_tract_weights.fig');
    close all;

end

end