% run_E4_real  E4: chb_hfov_porcine (HFOV 10 Hz, 44 frames/s), no ground truth.
% CHB_pig06-55.get = 60 s HFOV at 10 Hz (CHB_pig06-ref.get = conventional ventilation).
% Data: eidors3d.sourceforge.net/data_contrib/chb_hfov_porcine (cite Wolf et al., PCCM 2010).
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;

% E4a-E4c already run on 2026-11-16 (uncomment to repeat)
% cfg = kt_config();
% cfg.real.files = {fullfile(pwd, 'Data', 'chb_hfov_porcine', 'CHB_pig06-55.get')};
% cfg.real.t_range = [0 Inf];
% cfg.rec.timing_source = 'estimate';
% cfg.tag = 'E4a_hfov_full';
% cfg.notes = 'first look: whole 60 s HFOV record, MAP timing';
% kt_exp_real(cfg);

% E4b: is duty = 1.00 limited by the grid edge? (duty > 1 is not physical; a minimum above 1 flags a model problem)
% cfg = kt_config();
% cfg.real.files = {fullfile(pwd, 'Data', 'chb_hfov_porcine', 'CHB_pig06-55.get')};
% cfg.rec.timing_source = 'estimate';
% cfg.rec.duty_grid = 0:0.01:1.3;
% cfg.tag = 'E4b_dutygrid';
% cfg.notes = 'MAP duty hit 1.00 in E4a; extend grid to 1.3';
% kt_exp_real(cfg);

% E4c: 4th harmonic (41.3 Hz) folds to 2.72 Hz, near the pig heart rate; add an unaliased
% component at that bin and solve both jointly (real-data alias collision, cf. E3)
% cfg = kt_config();
% cfg.real.files = {fullfile(pwd, 'Data', 'chb_hfov_porcine', 'CHB_pig06-55.get')};
% cfg.rec.timing_source = 'estimate';
% cfg.real.extra_fold_of = 4;
% cfg.tag = 'E4c_h4_plus_cardiac';
% cfg.notes = 'kt_auto picked a non-physical alias for h4 in 11/12 blocks: test a 2.72 Hz line + h4 jointly';
% kt_exp_real(cfg);

% E4d-E4e already run on 2026-11-16 (uncomment to repeat)
% E4d: cardiac line at its own frequency (60 s spectrum: 2.67 Hz, h4 folds to 2.72 Hz;
% merged in 5 s blocks -> joint solve)
% cfg = kt_config();
% cfg.real.files = {fullfile(pwd, 'Data', 'chb_hfov_porcine', 'CHB_pig06-55.get')};
% cfg.rec.timing_source = 'estimate';
% cfg.real.cardiac_search = [1.5 3.5];
% cfg.tag = 'E4d_h4_cardiac_peak';
% cfg.notes = 'cardiac at detected peak instead of exactly at the h4 bin';
% kt_exp_real(cfg);

% E4e: duty optimum 1.19 > 1 in E4b: try voltages read one after another (per_meas)
% cfg = kt_config();
% cfg.real.files = {fullfile(pwd, 'Data', 'chb_hfov_porcine', 'CHB_pig06-55.get')};
% cfg.rec.timing_source = 'estimate';
% cfg.rec.timing_assumed.mode = 'per_meas';
% cfg.rec.duty_grid = 0:0.01:1.3;
% cfg.real.cardiac_search = [1.5 3.5];
% cfg.tag = 'E4e_per_meas';
% cfg.notes = 'per_meas schedule; does duty fall back to <= 1 and Phi drop?';
% kt_exp_real(cfg);


% E4f already run on 2026-11-16 (uncomment to repeat)
% E4f: one acquisition time per injection (equal-slot model gave duty > 1 in E4b/E4e)
% cfg = kt_config();
% cfg.real.files = {fullfile(pwd, 'Data', 'chb_hfov_porcine', 'CHB_pig06-55.get')};
% cfg.rec.timing_source = 'estimate';
% cfg.rec.timing_method = 'free';
% cfg.real.cardiac_search = [1.5 3.5];
% cfg.tag = 'E4f_free_timing';
% cfg.notes = 'free per-injection timing + cardiac line; block tol fixed (joint h4/cardiac in blocks)';
% kt_exp_real(cfg);


% E4g: held-out check of the schedule (estimate on one 30 s half, test on the other)
cfg = kt_config();
cfg.real.files = {fullfile(pwd, 'Data', 'chb_hfov_porcine', 'CHB_pig06-55.get')};
cfg.rec.timing_source = 'estimate';
cfg.rec.timing_method = 'free';
cfg.real.cardiac_search = [1.5 3.5];
cfg.real.cv_timing = true;
cfg.tag = 'E4g_timing_cv';
cfg.notes = 'two-fold CV of equal-slot and free timing';
kt_exp_real(cfg);
