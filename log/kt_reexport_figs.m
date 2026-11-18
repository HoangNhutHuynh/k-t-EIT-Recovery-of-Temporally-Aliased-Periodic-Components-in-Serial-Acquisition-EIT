function kt_reexport_figs(folder)
%KT_REEXPORT_FIGS Re-save every .fig under a folder as a white-background PNG.
%   kt_reexport_figs('runs')   % overwrites the .png next to each .fig
%   Use for runs saved before kt_light_fig existed (no need to rerun).
%
% Author: Hoang Nhut Huynh

f = dir(fullfile(folder, '**', '*.fig'));
for k = 1:numel(f)
    p = fullfile(f(k).folder, f(k).name);
    fh = openfig(p, 'invisible');
    kt_light_fig(fh);
    print(fh, [p(1:end-4) '.png'], '-dpng', '-r300');
    savefig(fh, p);
    close(fh);
    fprintf('re-exported %s\n', p);
end
end
