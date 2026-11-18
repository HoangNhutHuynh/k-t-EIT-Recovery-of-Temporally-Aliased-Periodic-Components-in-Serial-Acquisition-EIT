function run = kt_run_begin(cfg, experiment)
%KT_RUN_BEGIN Create a run folder, start the console diary and store the config.
%   run = kt_run_begin(cfg, experiment)
%   Folder: <cfg.out_root>/<yyyymmdd_HHMMSS>_<experiment>_<cfg.tag>/
%   Writes config.mat and config.txt (every field, one per line) there.
%   Always close with kt_run_end (also after an error) so the master logs
%   record failed attempts too.
%
% Author: Hoang Nhut Huynh

run.experiment = experiment;
run.cfg = cfg;
run.root = cfg.out_root;
run.t0 = tic;
base = sprintf('%s_%s_%s', datestr(now, 'yyyymmdd_HHMMSS'), experiment, kt_safe(cfg.tag));
run.id = base;
k = 1;
while exist(fullfile(run.root, run.id), 'dir')
    k = k + 1;
    run.id = sprintf('%s_%d', base, k);
end
run.dir = fullfile(run.root, run.id);
mkdir(run.dir);
run.fig_visible = get(0, 'DefaultFigureVisible');
if ~cfg.save_figs, set(0, 'DefaultFigureVisible', 'off'); end
diary off;
diary(fullfile(run.dir, 'console.txt'));
fprintf('=== k-t EIT run %s ===\n', run.id);
save(fullfile(run.dir, 'config.mat'), 'cfg');
lines = kt_flatten(cfg, '');
fid = fopen(fullfile(run.dir, 'config.txt'), 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);
fprintf('%s\n', lines{:});
end
