function kt_run_end(run, rows, summary, status)
%KT_RUN_END Write the metrics, append the master logs and close the diary.
%   kt_run_end(run, rows, summary, status)
%   rows    : n x 6 cell {param, pvalue, component, method, metric, value}
%   summary : one line with the key numbers of the run
%   status  : 'OK' (default) or 'FAILED: <message>'
%   Files written
%     <run.dir>/metrics.csv      metrics of this run (long format)
%     <root>/all_metrics.csv     the same rows for every run (kt_leaderboard)
%     <root>/run_log.csv         one line per run: id, status, key config, summary
%     <root>/RUN_LOG.md          a journal entry to complete by hand
%
% Author: Hoang Nhut Huynh

if nargin < 4, status = 'OK'; end
if isempty(rows), rows = cell(0, 6); end
elapsed = toc(run.t0);
cfg = run.cfg;
hdr = 'run_id,experiment,tag,param,pvalue,component,method,metric,value';
tag = kt_safe(cfg.tag);

lines = cell(size(rows, 1), 1);
for i = 1:size(rows, 1)
    lines{i} = sprintf('%s,%s,%s,%s,%.8g,%s,%s,%s,%.8g', run.id, run.experiment, tag, ...
        kt_safe(rows{i,1}), rows{i,2}, kt_safe(rows{i,3}), kt_safe(rows{i,4}), ...
        kt_safe(rows{i,5}), rows{i,6});
end
append_lines(fullfile(run.dir, 'metrics.csv'), hdr, lines);
append_lines(fullfile(run.root, 'all_metrics.csv'), hdr, lines);

keys = cfg.log_keys;
kv = cell(1, numel(keys));
for i = 1:numel(keys)
    kv{i} = sprintf('%s=%s', keys{i}, kt_val2str(kt_get(cfg, keys{i})));
end
key_cfg = strjoin(kv, '; ');
q = @(s) ['"' strrep(s, '"', '''') '"'];
append_lines(fullfile(run.root, 'run_log.csv'), ...
    'run_id,date,experiment,tag,status,elapsed_s,key_config,summary,notes', ...
    {sprintf('%s,%s,%s,%s,%s,%.1f,%s,%s,%s', run.id, datestr(now, 'yyyy-mm-dd HH:MM'), ...
    run.experiment, tag, q(status), elapsed, q(key_cfg), q(summary), q(cfg.notes))});

md = fullfile(run.root, 'RUN_LOG.md');
fid = fopen(md, 'a', 'n', 'UTF-8');
fprintf(fid, '\n## %s\n', run.id);
fprintf(fid, '- Status: %s (%.0f s)\n', status, elapsed);
fprintf(fid, '- Notes at launch (cfg.notes): %s\n', cfg.notes);
fprintf(fid, '- Key config: %s\n', key_cfg);
fprintf(fid, '- Key result: %s\n', summary);
fprintf(fid, '- Observations (fill in):\n- Decision / change for next run (fill in):\n');
fclose(fid);

fprintf('Run %s finished: %s (%.1f s)\nSummary: %s\nFolder: %s\n', ...
    run.id, status, elapsed, summary, run.dir);
set(0, 'DefaultFigureVisible', run.fig_visible);
diary off;
end

function append_lines(fname, hdr, lines)
newfile = ~exist(fname, 'file');
fid = fopen(fname, 'a');
if fid < 0, error('kt_run_end: cannot open %s', fname); end
if newfile, fprintf(fid, '%s\n', hdr); end
for i = 1:numel(lines)
    fprintf(fid, '%s\n', lines{i});
end
fclose(fid);
end
