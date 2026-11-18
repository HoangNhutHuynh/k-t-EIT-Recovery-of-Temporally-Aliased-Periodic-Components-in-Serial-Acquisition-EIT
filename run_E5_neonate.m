% run_E5_neonate  E5: second real record, ventilated neonate (if-neonate-ventilated, nm18).
% Goe-MF II at 13 frames/s, 208 stored measurements, 58 s. Heart rate ~2.47 Hz (0.19 cycles/frame):
% 2nd cardiac harmonic 4.9 Hz is below Nyquist (6.5 Hz), 3rd harmonic 7.4 Hz folds to 5.6 Hz (ell = -1).
% Ventilation ~0.49 Hz and its harmonics are unaliased and kept out of the noise floor.
% Weaker test than E4 (cardiac harmonics are small); its main value is a second subject for
% the schedule estimate and the alias selection.
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;

cfg = kt_config();
cfg.real.files = {fullfile(pwd, 'Data', 'if-neonate-ventilated', 'nm18.get')};
cfg.real.fs = 13;                       % .prl header; check console: f0 should be a plausible heart rate
cfg.real.f0_search = [2.45 2.6];         % cardiac fundamental (narrow: 5th ventilation harmonic is at 2.44 Hz)
cfg.real.n_harm = 3;
cfg.real.extra_excl = 0.4875*(1:12);     % ventilation and harmonics (0.0375 cycles/frame x 13)
cfg.rec.timing_source = 'estimate';
cfg.rec.timing_method = 'free';
cfg.real.cv_timing = true;
cfg.tag = 'E5_neonate_nm18';
cfg.notes = 'second subject: neonate cardiac harmonics, 13 fps';
kt_exp_real(cfg);
