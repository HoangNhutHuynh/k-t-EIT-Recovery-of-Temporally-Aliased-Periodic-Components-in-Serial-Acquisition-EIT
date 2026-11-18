function T = kt_leaderboard(root, metric, component, method, top)
%KT_LEADERBOARD Rank all logged runs by one metric (reads <root>/all_metrics.csv).
%   T = kt_leaderboard(root, metric, component, method, top)
%   metric    e.g. 'rho_ref', 'rho_true', 'dPhi_vs_yerworth'
%   component e.g. 'hfov_h3' ('' = all)
%   method    e.g. 'kt' ('' = all)
%   top       number of rows printed (default 15)
%   Metrics where smaller is better (relerr*, quad_tv*, phi_rel) are sorted
%   ascending, all others descending.
%   Example: kt_leaderboard('runs', 'rho_ref', 'hfov_h3', 'kt')
%
% Author: Hoang Nhut Huynh

if nargin < 3, component = ''; end
if nargin < 4, method = ''; end
if nargin < 5, top = 15; end
fid = fopen(fullfile(root, 'all_metrics.csv'), 'r');
if fid < 0, error('kt_leaderboard: %s/all_metrics.csv not found', root); end
c = textscan(fid, '%s %s %s %s %f %s %s %s %f', 'Delimiter', ',', 'HeaderLines', 1);
fclose(fid);
T.run_id = c{1}; T.tag = c{3}; T.param = c{4}; T.pvalue = c{5};
T.component = c{6}; T.method = c{7}; T.metric = c{8}; T.value = c{9};
keep = strcmp(T.metric, metric);
if ~isempty(component), keep = keep & strcmp(T.component, component); end
if ~isempty(method), keep = keep & strcmp(T.method, method); end
f = fieldnames(T);
for i = 1:numel(f), T.(f{i}) = T.(f{i})(keep); end
smaller_better = ~isempty(regexp(metric, '^(relerr|quad_tv|phi_rel)', 'once'));
if smaller_better
    [~, o] = sort(T.value, 'ascend');
else
    [~, o] = sort(T.value, 'descend');
end
for i = 1:numel(f), T.(f{i}) = T.(f{i})(o); end
fprintf('%-45s %-16s %-28s %-10s %-12s %10s\n', 'run_id', 'tag', 'param=pvalue', 'component', 'method', metric);
for i = 1:min(top, numel(T.value))
    fprintf('%-45s %-16s %-28s %-10s %-12s %10.4f\n', T.run_id{i}, T.tag{i}, ...
        sprintf('%s=%g', T.param{i}, T.pvalue(i)), T.component{i}, T.method{i}, T.value(i));
end
end
