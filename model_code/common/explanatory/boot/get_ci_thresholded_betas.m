function beta = get_ci_thresholded_betas(beta, ci)

% Threshold model betas based on zero-crossing sign differences between lower and upper CIs
% Joseph Griffis 2023

sign_l = sign(ci(1,:)); % lower CI signs
sign_u = sign(ci(2,:)); % upper CI signs
sign_dif = sign_l-sign_u; % difference in CI signs
beta(sign_dif~=0)=0; % remove betas with zero-crossing CIs
beta(beta == 0)=0; % NAN betas with zero values

end