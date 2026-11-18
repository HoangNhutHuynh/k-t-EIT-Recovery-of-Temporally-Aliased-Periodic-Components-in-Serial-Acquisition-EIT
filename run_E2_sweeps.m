% run_E2_sweeps  E2: robustness sweeps (linear forward model, one logged run per sweep).
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;

% S1: error in the assumed duty (true duty = 0.9) -> how exact must the timing be?
cfg = kt_config();
cfg.tag = 'S1';
cfg.rec.timing_source = 'assumed';
cfg.notes = 'timing mismatch, true duty 0.9';
kt_exp_sweep(cfg, 'rec.timing_assumed.duty', [0.5 0.7 0.8 0.85 0.9 0.95 1.0], 'rho_ref');

% S2: noise level
cfg = kt_config();
cfg.tag = 'S2';
cfg.notes = 'SNR sweep, oracle timing';
kt_exp_sweep(cfg, 'sim.snr_db', [20 30 40 50 60], 'rho_ref');

% S3: HFOV frequency (f0 = 11 Hz = fs/4 makes 33 Hz collide with -11 Hz; 22 Hz is Nyquist)
cfg = kt_config();
cfg.tag = 'S3';
cfg.notes = 'f0 sweep 6-20 Hz, oracle timing; rho_true is the valid score at collisions';
kt_exp_sweep(cfg, 'sim.f0', 6:1:20, 'rho_true');

% S4: regularisation (same alpha for every method)
cfg = kt_config();
cfg.tag = 'S4';
cfg.notes = 'alpha sweep';
kt_exp_sweep(cfg, 'rec.alpha', [0.001 0.003 0.01 0.03 0.1 0.3 1], 'rho_true');
