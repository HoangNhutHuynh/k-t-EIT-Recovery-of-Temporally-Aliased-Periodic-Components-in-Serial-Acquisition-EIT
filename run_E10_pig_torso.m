% run_E10_pig_torso  E10: six more Goe-MF II records (2 pigs, Data/Pig_torso), 13 -> 6.5 fps.
% Each record has ventilation at ~0.34 Hz (harmonics up to h9 stay below 3.25 Hz and only
% support the schedule estimate) and a strong cardiac line at 2.6-3.7 Hz. After removal of
% every second frame the cardiac line (g10-c005/c019/c029, g08-c035) or its 2nd harmonic
% (g08-c003/c017, added as 'line1') is aliased, while the full-rate image is the reference.
% Same settings as run_FINAL_E7 (pig_23kg mesh, time_prior baseline, schedule re-estimated
% from the reduced data). About 3-4 min per record. Metrics only, no figures.
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;
base = kt_config();
base.save_figs = false;
base.save_mat = false;
base.rec.eidors_temporal = {'time_prior'};
base.real.fs = 13;
base.real.f0_search = [0.30 0.42];      % ventilation fundamental
base.real.n_harm = 9;                   % h9 < 3.25 Hz: schedule support only
base.model.shape = 'pig_23kg';

% record, cardiac band (Hz), band of the 2nd cardiac harmonic (empty if above 6.5 Hz), n_harm
R = { 'g08-c003', [2.80 3.05], [5.75 5.95], 9;
      'g08-c017', [2.45 2.70], [5.05 5.25], 9;
      'g08-c035', [3.15 3.40], [],          8;    % f0 = 0.38 Hz here, h9 would be aliased
      'g10-c005', [3.40 3.65], [],          9;
      'g10-c019', [3.22 3.42], [],          9;
      'g10-c029', [3.60 3.76], [],          9 };

for r = 1:size(R, 1)
    cfg = base; cfg.tag = ['E10_' R{r, 1}];
    cfg.real.files = {fullfile(pwd, 'Data', 'Pig_torso', [R{r, 1} '.get'])};
    cfg.real.cardiac_search = R{r, 2};
    if ~isempty(R{r, 3}), cfg.real.extra_search = {R{r, 3}}; end
    cfg.real.n_harm = R{r, 4};
    cfg.notes = 'E10: Pig_torso cardiac line, 13 -> 6.5 fps';
    kt_exp_decim(cfg);
end
