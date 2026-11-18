% run_export_fig_data  Export the plotted data of the manuscript image figures (Figs. 3 to 5)
% to plain .mat files (*_data.mat next to each .fig), used by figures/fig_images.py.
% Reads the .fig files of the newest matching runs; nothing is reconstructed again.
% Requires run_FINAL_E7 and run_E10b_pig_torso_timing (figure section) to have been run.
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;
F = { '*_E7decim_FINAL_E7_p1130122', {'fig_decim_error_main', 'fig_decim_main'};
      '*_E7decim_FINAL_E7_hfov',     {'fig_decim_error', 'fig_decim_aligned'};
      '*_E7decim_FINAL_E7_neonate',  {'fig_decim_error', 'fig_decim_aligned'};
      '*_E7decim_E10_fig_g10-c005',  {'fig_decim_error', 'fig_decim_aligned'} };
for i = 1:size(F, 1)
    d = dir(fullfile(pwd, 'runs', F{i, 1}));
    d = d([d.isdir]);
    if isempty(d)
        fprintf('no run matching %s: skipped\n', F{i, 1});
        continue
    end
    [~, k] = max([d.datenum]);
    for n = F{i, 2}
        kt_export_fig_data(fullfile(d(k).folder, d(k).name, [n{1} '.fig']));
    end
end
