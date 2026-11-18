function kt_save_fig(run, fh, name)
%KT_SAVE_FIG Save a figure into the run folder as PNG (and .fig in MATLAB).
%   With cfg.save_figs = false the figure is only closed (nothing written).
%
% Author: Hoang Nhut Huynh

if ~run.cfg.save_figs
    close(fh);
    return
end
kt_light_fig(fh);
print(fh, fullfile(run.dir, [name '.png']), '-dpng', '-r300');
if ~exist('OCTAVE_VERSION', 'builtin')
    savefig(fh, fullfile(run.dir, [name '.fig']));
end
if ~run.cfg.keep_figs_open
    close(fh);
end
end
