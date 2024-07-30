function X = apply_dtlvc(X)

% Transform X to have unit norm of 1 
% i.e., apply direct total lesion volume control
% see Zhang et al., 2014 (Human Brain Mapping)
% Joseph Griffis 2023

% Divide each subject's lesion by the square root of the sum of lesion volume
X = X ./ sqrt(sum(X,2));

% Remove any NaNs/Infs
X(isnan(X))=0;
X(isinf(X))=0;
    
end