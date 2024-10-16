function cfg = check_and_fill_cfg_fields(cfg)
    
% This funciton checks for required fields in the cfg and fills them if missing
    
    % Model type
    if ~isfield(cfg, 'model_spec')
       disp('cfg.model_spec field is empty; setting to "plsr" with default settings')
       cfg.model_spec = 'plsr';
    end
    
    % Inferential modeling
    if ~isfield(cfg, 'fit_explanatory_model')
        disp('cfg.fit_explanatory_model field is empty; setting to 1 with default settings')
        cfg.fit_explanatory_model = 1;
    end

    % Hyper-parameter optimization
    if ~isfield(cfg, 'optimize_hyperparams') && ~contains(cfg.model_spec, {'municorr', 'bmunz', 'ttest', 'munilr', 'muniolsr', 'prop_sub'})
        disp('cfg.optimize_hyperparams field is empty; setting to 1 with default settings')
        cfg.optimize_hyperparams = 1;
    elseif ~isfield(cfg, 'optimize_hyperparams') && contains(cfg.model_spec, {'municorr', 'bmunz', 'ttest', 'munilr', 'muniolsr', 'prop_sub'})
        disp('Mass univariate analysis selected - hyper-parameter optimization will not be performed');
        cfg.optimize_hyperparams = 0;
    end
    
    % Cross-validation
    if ~isfield(cfg, 'cross_validation') && ~contains(cfg.model_spec, {'municorr', 'bmunz', 'ttest', 'munilr', 'muniolsr', 'prop_sub'})
        disp('cfg.cross_validation field is empty; setting to 1 with default settings')
        cfg.cross_validation = 1;
    elseif ~isfield(cfg, 'cross_validation') && contains(cfg.model_spec, {'municorr', 'bmunz', 'ttest', 'munilr', 'muniolsr', 'prop_sub'})
        disp('Mass univariate analysis selected - cross-validation will not be performed');
        cfg.cross_validation = 0;
    end
    
    % Bootstrapping
    if ~isfield(cfg, 'bootstrap') && ~contains(cfg.model_spec, {'municorr', 'bmunz', 'ttest', 'munilr', 'muniolsr', 'prop_sub'})
        disp('cfg.bootstrap field is empty; setting to 1 with default settings')
        cfg.bootstrap = 1;
    elseif ~isfield(cfg, 'bootstrap') && contains(cfg.model_spec, {'municorr', 'bmunz', 'ttest', 'munilr', 'muniolsr', 'prop_sub'})
        disp('Mass univariate analysis selected - bootstrapping is disabled')
        cfg.bootstrap = 0;    
    end
    
    % Permutation testing
    if ~isfield(cfg, 'permutation')
        disp('cfg.permutation field is empty; setting to 1 with default settings')
        cfg.permutation = 1;
    end
    
    % Jack-knife testing
    if ~isfield(cfg, 'jackknife')
         disp('cfg.jackknife field is empty; setting to 0 with default settings')
         cfg.jackknife = 0;   
    end

end