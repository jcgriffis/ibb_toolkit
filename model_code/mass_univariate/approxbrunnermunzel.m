function [o]=approxbrunnermunzel(x, y, sidedness)
% sidedness = 1 -> two sided; = 2 -> greater ; = 3 -> less
% ported from
% https://github.com/toshi-ara/brunnermunzel/blob/master/src/bm_test.f

alpha = 0.05; % for computing 95% CI
TOLER = 10^(-10);

P1 = [0.0, 1.0, 0.0];
P0 = [0.0, 0.0, 1.0];

% columnize data
x = x(:);
y = y(:);
% remove nan if any
nonnanidx = find(~isnan(x));
x = x(nonnanidx); clear nonnanidx
nonnanidx = find(~isnan(y));
y = y(nonnanidx); clear nonnanidx

xy = [x;y];

nx = length(x);
ny = length(y);
nxy = nx+ny;
[rkx, tiedadjx] = tiedrank(x);
[rky, tiedadjy] = tiedrank(y);
[rkxy, tiedadjxy] = tiedrank(xy);

mx = mean(rkxy(1:nx));
my = mean(rkxy(nx+1:nxy));

pst = (my - (ny + 1)*0.5)/nx; % pst: estimation of "P(X<Y)+.5*P(X=Y)"

if abs(pst-1) < TOLER       % pst == 1, non-overlapped data: X < Y
    ci(1:2) = pst;          % (/1.0, 1.0/)
    stat_ = Inf;            % Inf
    df = NaN;               % NaN
    pval = P1(sidedness);   % P = 1 in "greater", P = 0 in others
    se = NaN; ci=[NaN NaN];
elseif abs(pst) < TOLER     % pst == 0, non-overlapped data: X > Y
    ci(1:2) = pst;          % (/0.0, 0.0/)
    stat_ = -Inf;           % -Inf
    df = NaN;               % NaN
    pval = P0(sidedness);   % P = 1 in "less", P = 0 in others
    se = NaN; ci=[NaN NaN];
else                        % overlapped data
    [stat_, df, se] = calc_stat(nx, ny, rkx, rky, rkxy, mx, my);
    pval = calc_pval(stat_, df, sidedness);
    ci = calc_confint(pst, df, se, alpha);
end
o.pval=pval; o.stat=stat_;o.df=df;o.pst=pst;o.se=se;o.ci=ci;
% text output
fprintf(1,'\n Brunner-Munzel Test Statistic = %f, df = %f, p-value = %f',stat_, df, pval)
fprintf(1,'\n sample estimates of P(X<Y)+.5*P(X=Y): %f', pst)
fprintf(1,'\n 95 percent confidence interval: %f, %f',ci)
fprintf(1,'\n');
return;

function [stat_, df, se] = calc_stat(nx, ny, rkx, rky, rkxy, mx, my)

n1 = nx; n2 = ny;
vx = 0; vy = 0;
dx = (rkxy(1:nx) - rkx - mx + (nx + 1) * 0.5).^2;
dy = (rkxy(nx+1:nx+ny) - rky - my + (ny + 1) * 0.5).^2;
vx = sum(dx(1:nx));
vy = sum(dy(1:ny));

% variance of group x and y
vx = vx ./(nx - 1);
vy = vy ./(ny - 1);

nvx = n1 .* vx;
nvy = n2 .* vy;
nv = nvx + nvy;

stat_ = n1 * n2 ./(nx + ny) .*(my - mx) ./ sqrt(nv);
df = nv .* nv ./(nvx .* nvx ./(nx - 1) + nvy .* nvy ./(ny - 1));
se = sqrt(vx ./(n1 .* n2 .* n2) + vy ./(n1 .* n1 .* n2));
return


function pval = calc_pval(stat_, df, sidedness)
lowertail = [0, 1, 0];
multi = [2.0, 1.0, 1.0];


if sidedness == 1
    pval = tcdf(abs(stat_), df, 'upper') * multi(sidedness);
elseif sidedness == 2
    pval = (tcdf(stat_, df)) * multi(sidedness);
elseif sidedness == 3
    pval = (1-tcdf(stat_, df)) * multi(sidedness);
else
    fprintf(1,'\nERROR in calc_pval, return NaN\n');
    pval = NaN;
end

return;

function ci = calc_confint(pst, df, se, alpha)

ci = pst + tinv(alpha * 0.5, abs(df)) * se * [-1 1] ;
ci(1) = max([0, ci(1)]);
ci(2) = min([1, ci(2)]);
ci = sort(ci);

return





