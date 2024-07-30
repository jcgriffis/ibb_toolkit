function opt_thresh = get_optimal_threshold(obs_y, roc_obj, cost)

idx_1 = roc_obj.Metrics.ClassName == 1;
X = roc_obj.Metrics(idx_1,:).FalsePositiveRate;
Y = roc_obj.Metrics(idx_1,:).TruePositiveRate;
T = roc_obj.Metrics(idx_1,:).Threshold;

p = sum(obs_y == 1);
n = sum(obs_y == -1);
m = (cost(2,1)-cost(2,2))/(cost(1,2)-cost(1,1))*n/p;
[~,idx_opt] = min(X - Y/m);
opt_thresh = T(idx_opt);

end