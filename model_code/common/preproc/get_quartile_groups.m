function [groups] = get_quartile_groups(strat_var)
% Stratify data into 4 groups based on quartiles of a stratification variable
% Joseph Griffis 2023

% Get quartile boundaries
q1thr = prctile(strat_var, 25);
q2thr = prctile(strat_var, 50);
q3thr = prctile(strat_var, 75);

% Find patients within each quartile
s1 = find(strat_var <= q1thr);
s2 = find(strat_var > q1thr & strat_var <= q2thr);
s3 = find(strat_var > q2thr & strat_var <= q3thr);
s4 = find(strat_var > q3thr);
disp(['group 1:' num2str(numel(s1))]);
disp(['group 2:' num2str(numel(s2))]);
disp(['group 3:' num2str(numel(s3))]);
disp(['group 4:' num2str(numel(s4))]);

% Generate grouping variable
groups = zeros(length(strat_var),1);
groups(s1) = 1;
groups(s2) = 2;
groups(s3) = 3;
groups(s4) = 4;

end