function write_model_weights_connectome(model_results)

%%% output .node and .edge files for external viewers (e.g. MRIcroGL)

% navigate to output directory
cd(model_results.cfg.out_dir);

if ~strcmp(model_results.cfg.model_spec, 'municorr')
    write_connectome_weights_mlsm(model_results);
else
    write_connectome_weights_ulsm(model_results);
end

end