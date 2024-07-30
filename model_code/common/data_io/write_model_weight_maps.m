function write_model_weight_maps(model_results)

% Writes out frequency map and model weight maps for voxel-based lesion analyses; probably needs to be refactored because there's a lot going on
% Joseph Griffis 2024

% Go to output directory
cd(model_results.cfg.out_dir);

% Get mask image and header info
nifti_out = niftiread(model_results.cfg.mask_path);
nifti_hdr = niftiinfo(model_results.cfg.mask_path);

% Change header info for continuous output images
nifti_hdr.Datatype = 'single';
nifti_hdr.BitsPerPixel = 32;
nifti_hdr.raw.dime.datatype= 16;
nifti_hdr.raw.dime.bitpix= 32;
nifti_out = single(nifti_out);

% Mass univariate correlations
if strcmp(model_results.cfg.model_spec, 'municorr') || strcmp(model_results.cfg.model_spec, 'munilr') || strcmp(model_results.cfg.model_spec, 'bmunz') || strcmp(model_results.cfg.model_spec, 'ttest')
    write_mass_unicorr_weight_maps(nifti_out, nifti_hdr, model_results);
else
    write_mlsm_weight_maps(nifti_out, nifti_hdr, model_results);
end

end