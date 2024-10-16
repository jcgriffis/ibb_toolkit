function write_model_weights_connectome(model_results)

%%% output .node and .edge files for external viewers (e.g. MRIcroGL)

% navigate to output directory
cd(model_results.cfg.out_dir);

if strcmp(model_results.cfg.model_spec, 'municorr') || strcmp(model_results.cfg.model_spec, 'munilr') || strcmp(model_results.cfg.model_spec, 'bmunz') || strcmp(model_results.cfg.model_spec, 'ttest') || strcmp(model_results.cfg.model_spec, 'muniolsr') || strcmp(model_results.cfg.model_spec, 'prop_sub')
    write_connectome_weights_ulsm(model_results);
else
    write_connectome_weights_mlsm(model_results);
end

end