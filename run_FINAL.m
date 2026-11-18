% run_FINAL  Final runs for the manuscript, with figures and results.mat (about 40-60 min).
% Run only after E7 (decimation) and E8 (thorax mesh + baselines) have been checked.
% Configuration frozen on 2026-11-16 after E1-E6: MAP timing (equal slots) for simulation,
% free per-injection timing + two-fold CV for real data, cardiac line joint with h4 for HFOV.
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;
base = kt_config();
base.save_figs = true;
base.save_mat = true;

% --- simulation
cfg = base; cfg.tag = 'FINAL_E1a_oracle'; cfg.notes = 'final: oracle timing';
kt_exp_sim(cfg);

cfg = base; cfg.tag = 'FINAL_E1d_lag90'; cfg.rec.timing_source = 'estimate';
for c = 1:numel(cfg.sim.components)
    h = cfg.sim.components(c).harm;
    if h > 0, cfg.sim.components(c).phase_lr_deg = 90*h; end
end
cfg.notes = 'final: 90 deg left-right lag, MAP timing';
kt_exp_sim(cfg);

cfg = base; cfg.tag = 'FINAL_S1'; cfg.rec.timing_source = 'assumed'; cfg.notes = 'final: timing mismatch';
kt_exp_sweep(cfg, 'rec.timing_assumed.duty', [0.5 0.7 0.8 0.85 0.9 0.95 1.0], 'rho_ref');

cfg = base; cfg.tag = 'FINAL_S3'; cfg.notes = 'final: f0 sweep';
kt_exp_sweep(cfg, 'sim.f0', 6:1:20, 'rho_true');

cfg = base; cfg.tag = 'FINAL_E3'; cfg.notes = 'final: alias collision f0 = fs/4';
kt_exp_collision(cfg, [0.25 0.5 0.75 0.9 1.0]);

% --- real data
rbase = base;
rbase.rec.timing_source = 'estimate';
rbase.rec.timing_method = 'free';
rbase.real.cv_timing = true;
rbase.rec.eidors_temporal = {'time_prior'};

cfg = rbase; cfg.tag = 'FINAL_E4_hfov';
cfg.real.files = {fullfile(pwd, 'Data', 'chb_hfov_porcine', 'CHB_pig06-55.get')};
cfg.real.cardiac_search = [1.5 3.5];
cfg.model.shape = 'pig_23kg';
cfg.notes = 'final: HFOV pig';
kt_exp_real(cfg);

cfg = rbase; cfg.tag = 'FINAL_E5_neonate';
cfg.real.files = {fullfile(pwd, 'Data', 'if-neonate-ventilated', 'nm18.get')};
cfg.real.fs = 13; cfg.real.f0_search = [2.45 2.6]; cfg.real.n_harm = 3;
cfg.real.extra_excl = 0.4875*(1:12);
cfg.model.shape = 'neonate';
cfg.notes = 'final: neonate';
kt_exp_real(cfg);

files = {'p1130107.get', 'p1130122.get'};
for i = 1:numel(files)
    cfg = rbase; cfg.tag = ['FINAL_E6_' strtok(files{i}, '.')];
    cfg.real.files = {fullfile(pwd, 'Data', 'pig_2003', files{i})};
    cfg.real.fs = 13; cfg.real.f0_search = [0.6 0.8]; cfg.real.n_harm = 12; cfg.real.block_sec = 10;
    cfg.model.shape = 'pig_23kg';
    cfg.notes = 'final: pig_2003';
    kt_exp_real(cfg);
end

% --- frame-rate reduction (real-data ground truth)
for f = {'p1130107', 'p1130122'}
    cfg = base; cfg.tag = ['FINAL_E7_' f{1}];
    cfg.rec.eidors_temporal = {'time_prior'};
    cfg.real.files = {fullfile(pwd, 'Data', 'pig_2003', [f{1} '.get'])};
    cfg.real.fs = 13; cfg.real.f0_search = [0.6 0.8]; cfg.real.n_harm = 9;
    cfg.model.shape = 'pig_23kg';
    cfg.notes = 'final: pig_2003 13 -> 6.5 fps';
    kt_exp_decim(cfg);
end
cfg = base; cfg.tag = 'FINAL_E7_hfov';
cfg.rec.eidors_temporal = {'time_prior'};
cfg.real.files = {fullfile(pwd, 'Data', 'chb_hfov_porcine', 'CHB_pig06-55.get')};
cfg.real.n_harm = 2; cfg.real.cardiac_search = [1.5 3.5];
cfg.model.shape = 'pig_23kg';
cfg.notes = 'final: HFOV 44 -> 22 fps';
kt_exp_decim(cfg);
cfg = base; cfg.tag = 'FINAL_E7_neonate';
cfg.rec.eidors_temporal = {'time_prior'};
cfg.real.files = {fullfile(pwd, 'Data', 'if-neonate-ventilated', 'nm18.get')};
cfg.real.fs = 13; cfg.real.f0_search = [2.45 2.6]; cfg.real.n_harm = 2;
cfg.real.extra_excl = 0.4875*(1:12);
cfg.model.shape = 'neonate';
cfg.notes = 'final: neonate 13 -> 6.5 fps';
kt_exp_decim(cfg);
