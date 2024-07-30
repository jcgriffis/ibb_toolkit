function cv_results = summarize_cv_results(cv_results, cfg, perm_flag)

if isfield(cfg.cv, 'summary_type')
    if strcmp(cfg.cv.summary_type, 'mean')
        if cfg.cat_Y == 0 && perm_flag == 0
            cv_results.avg.explained(:,1) = nanmean(cv_results.explained,2);
            cv_results.avg.r2_ss(:,1) = nanmean(cv_results.r2_ss,2);
            cv_results.avg.corr(:,1) = nanmean(cv_results.corr,2);
            cv_results.avg.mse(:,1) = nanmean(cv_results.mse, 2);
        elseif cfg.cat_Y == 1 && perm_flag == 0
            cv_results.avg.classrate(:,1) = nanmean(cv_results.classrate(:,:,1), 2);
            cv_results.avg.classrate(:,2) = nanmean(cv_results.classrate(:,:,2), 2);
            cv_results.avg.classrate(:,3) = nanmean(cv_results.classrate(:,:,3), 2);
            cv_results.avg.roc_auc(:,1) = nanmean(cv_results.roc_auc, 2);
        elseif cfg.cat_Y == 0 && perm_flag == 1
            cv_results.avg.mse = nanmean(cv_results.mse,2);
        elseif cfg.cat_Y == 1 && perm_flag == 1
            cv_results.avg.classrate(:,1) = nanmean(cv_results.classrate(:,:,1), 2);
            cv_results.avg.classrate(:,2) = nanmean(cv_results.classrate(:,:,2), 2);
            cv_results.avg.roc_auc(:,1) = nanmean(cv_results.roc_auc, 2);
        end
    elseif strcmp(cfg.cv.summary_type, 'median')
        if cfg.cat_Y == 0 && perm_flag == 0
            cv_results.avg.explained(:,1) = nanmedian(cv_results.explained,2);
            cv_results.avg.r2_ss(:,1) = nanmedian(cv_results.r2_ss,2);
            cv_results.avg.corr(:,1) = nanmedian(cv_results.corr,2);
            cv_results.avg.mse(:,1) = nanmedian(cv_results.mse, 2);
        elseif cfg.cat_Y == 1 && perm_flag == 0
            cv_results.avg.classrate(:,1) = nanmedian(cv_results.classrate(:,:,1), 2);
            cv_results.avg.classrate(:,2) = nanmedian(cv_results.classrate(:,:,2), 2);
            cv_results.avg.classrate(:,3) = nanmedian(cv_results.classrate(:,:,3), 2);
            cv_results.avg.roc_auc(:,1) = nanmedian(cv_results.roc_auc, 2);
        elseif cfg.cat_Y == 0 && perm_flag == 1
            cv_results.avg.mse = nanmedian(cv_results.mse,2);
        elseif cfg.cat_Y == 1 && perm_flag == 1
            cv_results.avg.classrate(:,1) = nanmedian(cv_results.classrate(:,:,1), 2);
            cv_results.avg.classrate(:,2) = nanmedian(cv_results.classrate(:,:,2), 2);
            cv_results.avg.roc_auc(:,1) = nanmedian(cv_results.roc_auc, 2);
        end
end

end