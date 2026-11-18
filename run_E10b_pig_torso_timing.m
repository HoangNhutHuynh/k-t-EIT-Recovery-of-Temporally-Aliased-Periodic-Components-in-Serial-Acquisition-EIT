% run_E10b_pig_torso_timing  E10b: schedule estimation with cross-validation on the Pig_torso
% records at the native 13 frames/s (same analysis as run_E6 / FINAL E4-E6, for Table II of the
% manuscript). Components: ventilation harmonics h1-h9 (0.35 Hz series) plus the cardiac line.
% One record per piglet is enough for the table; the other records are optional.
% Runtime: free timing + CV with 10 components, about 5-10 min per record.
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;

R = { 'g08-c017', [2.45 2.70];
      'g10-c005', [3.40 3.65] };
for r = 1:size(R, 1)
    cfg = kt_config();
    cfg.real.files = {fullfile(pwd, 'Data', 'Pig_torso', [R{r, 1} '.get'])};
    cfg.real.fs = 13;
    cfg.real.f0_search = [0.30 0.42];
    cfg.real.n_harm = 9;
    cfg.real.cardiac_search = R{r, 2};
    cfg.real.block_sec = 10;
    cfg.model.shape = 'pig_23kg';
    cfg.rec.timing_source = 'estimate';
    cfg.rec.timing_method = 'free';
    cfg.real.cv_timing = true;
    cfg.rec.eidors_temporal = {};
    cfg.tag = ['E10b_' R{r, 1}];
    cfg.notes = 'E10b: Pig_torso schedule estimation + CV at 13 fps';
    kt_exp_real(cfg);
end

% --- figure for the manuscript: error images of the aliased cardiac line, piglet 2 (g10-c005)
cfg = kt_config();
cfg.save_figs = true;
cfg.save_mat = false;
cfg.rec.eidors_temporal = {'time_prior'};
cfg.real.files = {fullfile(pwd, 'Data', 'Pig_torso', 'g10-c005.get')};
cfg.real.fs = 13; cfg.real.f0_search = [0.30 0.42]; cfg.real.n_harm = 9;
cfg.real.cardiac_search = [3.40 3.65];
cfg.model.shape = 'pig_23kg';
cfg.decim.fig_rows = {'cardiac'};
cfg.tag = 'E10_fig_g10-c005';
cfg.notes = 'E10: figure row, piglet 2 cardiac line 13 -> 6.5 fps';
kt_exp_decim(cfg);
