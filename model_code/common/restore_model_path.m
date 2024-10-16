function restore_model_path(model_dir)

% Reset model path to avoid conflicts and having too many unnecessary directories on path
% Joseph Griffis 2023

my_path = string(path);

if contains(my_path, fullfile(model_dir, 'reglm'))
    rmpath(genpath(fullfile(model_dir, 'reglm')));
end
if contains(my_path, fullfile(model_dir, 'nlinsvr'))
    rmpath(genpath(fullfile(model_dir, 'nlinsvr')));
end
if contains(my_path, fullfile(model_dir, 'plsr'))
    rmpath(genpath(fullfile(model_dir, 'plsr')));
end
if contains(my_path, fullfile(model_dir, 'mass_univariate'))
    rmpath(genpath(fullfile(model_dir, 'mass_univariate')));
end

end