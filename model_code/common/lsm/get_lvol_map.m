function get_lvol_map(X, X_inc_ind, cfg)

% Generate map showing correlations between lesion location and lesion volume
% Joseph Griffis 2024

if ~isfield(cfg, 'lvol')
    cfg.lvol = sum(X,2);
end

% Get bivariate correlation between lesion location and lesion volume
lvol_map = corr(X, cfg.lvol);

% Read in NIFTI file and header
nifti_out = niftiread(cfg.mask_path);
nifti_hdr = niftiinfo(cfg.mask_path);

% Change header info for continuous output images
nifti_hdr.Datatype = 'single';
nifti_hdr.BitsPerPixel = 32;
nifti_hdr.raw.dime.datatype= 16;
nifti_hdr.raw.dime.bitpix= 32;
nifti_out = single(nifti_out);

% Write out map
nifti_out(:) = 0;
nifti_out(cfg.mask_inds(X_inc_ind)) = lvol_map;
niftiwrite(nifti_out, 'lesion_volume_corr_map', nifti_hdr); 

end