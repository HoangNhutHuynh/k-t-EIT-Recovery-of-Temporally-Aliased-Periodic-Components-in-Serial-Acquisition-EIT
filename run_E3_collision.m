% run_E3_collision  E3: two aliases in one bin (f0 = fs/4), joint solve and separability.
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;

cfg = kt_config();
cfg.tag = 'E3_per_stim';
cfg.notes = 'f0 = 11 Hz, 33 Hz folds onto -11 Hz';
kt_exp_collision(cfg, [0.25 0.5 0.75 0.9 1.0]);

% same with voltages read one after another (more distinct acquisition times)
cfg = kt_config();
cfg.tag = 'E3_per_meas';
cfg.acq.timing = 'per_meas';
cfg.notes = 'per_meas schedule';
kt_exp_collision(cfg, [0.25 0.5 0.75 0.9 1.0]);
