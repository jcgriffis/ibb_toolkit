% This function projects voxel-level weight vectors back to the brain
% and applies thresholding and rescaling if desired.
% Joseph Griffis, 2017, Washington University in St Louis
function [temp] = get_lmwv(temp_nii, lm_has_dmg, lm_weights, rescale_opt, thresh, direction)

% Inputs
% temp_nii - nifti file with same dimensions and in same space as lesion data
% lm_has_damge - indices of brain mask corresponding with temp_nii for voxels used in analyses
% lm_weights - voxel weights for included voxels (same length as lm_has_damage)
% rescale_flag - flag to indicate rescaling (see below)
% thresh - threshold to apply if rescale_flag == 1
% direction - effect direction ('pos' or 'neg') if rescale_flag == 1

% Outputs
% nifti file containing rescaled and thresholded (if indicated) model weights

% Put back in NIFTI file
temp = temp_nii;
temp.img = single(temp.img);
temp.img(:,:,:)=0;
temp.img(lm_has_dmg) = lm_weights;

% Rescale as indicated
if ~isempty(rescale_opt)
    if strcmp(rescale_opt, 'percentile')
        if strcmp(direction, 'pos')==1
            temp.img(temp.img < 0) = 0;
            p_thresh = prctile(temp.img(temp.img > 0), thresh);
            temp.img(temp.img < p_thresh) = 0;
            temp.img(lm_has_dmg) = temp.img(lm_has_dmg)./max(temp.img(lm_has_dmg));
            temp.img = temp.img.*10;
        elseif strcmp(direction, 'neg')==1
            temp.img = temp.img.*-1;
            temp.img(temp.img < 0) = 0;
            p_thresh = prctile(temp.img(temp.img > 0), thresh);
            temp.img(temp.img < p_thresh) = 0;
            temp.img(lm_has_dmg) = temp.img(lm_has_dmg)./max(temp.img(lm_has_dmg));
            temp.img = temp.img.*10;
        else
            p_thresh = prctile(temp.img(temp.img ~= 0), thresh);
            temp.img(temp.img < p_thresh) = 0;
            temp.img(lm_has_dmg) = temp.img(lm_has_dmg)./max(temp.img(lm_has_dmg));
            temp.img = temp.img.*10;        
        end
    
    elseif strcmp(rescale_opt, 'zscore')
        temp.img(lm_has_dmg) = zscore(temp.img(lm_has_dmg));
    
    elseif strcmp(rescale_opt, 'max')
        if strcmp(direction, 'pos')
             temp.img(lm_has_dmg) = temp.img(lm_has_dmg)./max(temp.img(lm_has_dmg));
             temp.img = temp.img.*10;  
        elseif strcmp(direction, 'neg')
             temp.img(lm_has_dmg) = temp.img(lm_has_dmg)./min(temp.img(lm_has_dmg));
             temp.img = temp.img.*10;              
        else
             temp.img(lm_has_dmg) = temp.img(lm_has_dmg)./max(abs(temp.img(lm_has_dmg)));
             temp.img = temp.img.*10;            
        end
    elseif strcmp(rescale_opt, 'pval')
        temp.img(lm_has_dmg) = 1-temp.img(lm_has_dmg);
    elseif isempty(rescale_opt)
        if strcmp(direction, 'pos')
            temp.img(temp.img < 0) = 0;
        elseif strcmp(direction, 'neg')
            temp.img(temp.img > 0) = 0;
        end
    end
end

% Set header
temp.hdr.dime.datatype=16;
temp.hdr.dime.bitpix=32;
temp.hdr.dime.glmax = max(temp.img(:));
temp.hdr.dime.glmin = min(temp.img(:));

end