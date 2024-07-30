function model_results = get_boot_statistics(X, Y, stat, bootstat, bootfun, model_results, cfg)


if cfg.cat_Y == 1
    % Get bootstrapped p-values
    if cfg.boot.get_pvals == 1
        model_results = get_boot_pvals_cl(stat, bootstat, model_results, cfg);   
    end
    
    % Get CI for model classification rates
    if cfg.boot.get_cis == 1
        model_results = get_boot_cis_cl(stat, bootstat, bootfun, X ,Y, model_results, cfg);
    end
elseif cfg.cat_Y == 0
    % Get bootstrapped p-values
    if cfg.boot.get_pvals == 1
        model_results = get_boot_pvals_reg(stat, bootstat, model_results, cfg);   
    end
    
    % Get CI for model classification rates
    if cfg.boot.get_cis == 1
        model_results = get_boot_cis_reg(stat, bootstat, bootfun, X ,Y, model_results, cfg);
    end
end

% Save bootstrap statistics if indicated
if cfg.boot.save_boot_results == 1
    if cfg.parallel == 1
        bootstat = gather(bootstat);
        bootstat = bootstat(:,1);
    end
    if isfile('boot_results.mat')
       save('boot_results.mat', 'bootstat', '-append');
    else
       save('boot_results.mat', 'bootstat', '-v7.3');
    end
end

end