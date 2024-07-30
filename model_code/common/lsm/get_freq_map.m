function get_freq_map(X, X_inc_ind, cfg)

% Generates a frequency map for the predictor data
% Joseph Griffis 2024

freq_map = sum(X >= cfg.min_obs,1);
nifti_out = niftiread(cfg.mask_path);
nifti_hdr = niftiinfo(cfg.mask_path);
nifti_out(:) = 0;
nifti_out(cfg.mask_inds(X_inc_ind)) = freq_map;
niftiwrite(nifti_out, 'lesion_frequency_map', nifti_hdr); 

end