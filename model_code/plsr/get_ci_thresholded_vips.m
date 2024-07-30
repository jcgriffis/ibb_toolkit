function vip_thresh_l = get_ci_thresholded_vips(vip_scores, ci)
% Threshold PLSR VIPs based on whether the lower value of the CI crosses a
% threshold
% Joseph Griffis 2023

    vip_thresh_l = vip_scores; % copy VIP vector for thresholding
    my_cis_l = ci(1,:); % lower CI
    mean_vip = mean(vip_scores(:));
    vip_thresh_l(my_cis_l < mean_vip)=0; % remove betas with zero-crossing CIs
 
end