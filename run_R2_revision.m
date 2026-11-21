% run_R2_revision  Analyses for the strict pre-submission review (R2) of the k-t EIT manuscript.
% Sections can be run separately (Ctrl+Enter on a %% cell) or all at once.
% Metrics only, no figures. Every section appends its rows to one CSV in runs/:
%   runs/R2_<section>.csv  columns: tag, record, param, pvalue, component, method, metric, value
% Approximate run times (desktop PC): Y1 40 min, Y2 50 min, Y3 30 min, Y4 70 min,
% Y5 70 min, Y6 5 min, Y7 30 min, Y8 10 min.
%
%   Y1 (M1)  every baseline also reconstructed at the full rate: each method scored
%            against its own full-rate image (rho_self) and against the full-rate image
%            of every other method (rho_ref_<method>), nine records, alpha = 0.03
%   Y2 (M7)  schedule ablation for k-t and phase correction: FTC schedule ('uniform')
%            and the full-rate schedule ('full') instead of the re-estimated one
%   Y3 (M4)  unaliased nuisance line at the folded bin of every aliased component,
%            solved jointly by k-t (real-data collision with a measured reference)
%   Y4 (M8)  alpha = 0.01 and 0.1 with all four baselines, nine records
%   Y5 (M8)  Laplace and Tikhonov priors with all four baselines, nine records
%   Y6 (m3)  HFOV with the re-estimated duty bounded by 1/D
%   Y7 (m1)  stagger of the voltages inside an injection (kappa = 0 ... 1), native rate
%   Y8 (m12) simulation with model misspecification: stagger, timing jitter, f0 drift
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;
if ~exist('RUN', 'var'), RUN = {'Y1', 'Y2', 'Y3', 'Y4', 'Y5', 'Y6', 'Y7', 'Y8'}; end   % e.g. RUN = {'Y3'}; run_R2_revision

base = kt_config();
base.save_figs = false; base.save_mat = false;
base.rec.eidors_temporal = {'time_prior'};
base.real.extra_search = {};
R = kt_r2_records(base);
if ~exist('runs', 'dir'), mkdir('runs'); end   % fresh clone: runs/ is not tracked
open_csv = @(s) fopen(fullfile('runs', ['R2_' s '.csv']), 'a');

%% Y1: alternative full-rate references (M1)
if any(strcmp(RUN, 'Y1'))
    fid = open_csv('Y1_altref');
    for r = 1:numel(R)
        cfg = R(r).cfg; cfg.decim.alt_refs = true;
        cfg.tag = ['R2_Y1_' R(r).name]; cfg.notes = 'R2 Y1: alternative full-rate references';
        out = kt_exp_decim(cfg); kt_write_rows(fid, cfg.tag, R(r).name, out.rows);
    end
    fclose(fid);
end

%% Y2: schedule ablation (M7)
if any(strcmp(RUN, 'Y2'))
    fid = open_csv('Y2_schedule');
    for tm = {'uniform', 'full'}
        for r = 1:numel(R)
            cfg = R(r).cfg; cfg.decim.timing = tm{1};
            cfg.rec.methods = {'naive', 'yerworth', 'kt'}; cfg.rec.try_eidors_ftc = false; cfg.rec.eidors_temporal = {};
            cfg.tag = ['R2_Y2_' tm{1} '_' R(r).name]; cfg.notes = 'R2 Y2: schedule ablation';
            out = kt_exp_decim(cfg); kt_write_rows(fid, cfg.tag, R(r).name, out.rows);
        end
    end
    fclose(fid);
end

%% Y3: joint solve with an unaliased nuisance line at the folded bin (M4)
if any(strcmp(RUN, 'Y3'))
    fid = open_csv('Y3_nuisance');
    for r = 1:numel(R)
        cfg = R(r).cfg; cfg.decim.nuisance_l0 = true;
        cfg.rec.methods = {'yerworth', 'kt_single', 'kt'}; cfg.rec.try_eidors_ftc = false; cfg.rec.eidors_temporal = {};
        cfg.tag = ['R2_Y3_' R(r).name]; cfg.notes = 'R2 Y3: nuisance line at the folded bin, joint solve';
        out = kt_exp_decim(cfg); kt_write_rows(fid, cfg.tag, R(r).name, out.rows);
    end
    fclose(fid);
end

%% Y4: regularization level, all baselines (M8)
if any(strcmp(RUN, 'Y4'))
    fid = open_csv('Y4_alpha');
    for a = [0.01 0.1]
        for r = 1:numel(R)
            cfg = R(r).cfg; cfg.rec.alpha = a;
            cfg.tag = sprintf('R2_Y4_alpha%g_%s', a, R(r).name); cfg.notes = 'R2 Y4: alpha sweep, all baselines';
            out = kt_exp_decim(cfg); kt_write_rows(fid, cfg.tag, R(r).name, out.rows);
        end
    end
    fclose(fid);
end

%% Y5: other priors, all baselines (M8)
if any(strcmp(RUN, 'Y5'))
    fid = open_csv('Y5_prior');
    for pr = {'laplace', 'tikhonov'}
        for r = 1:numel(R)
            cfg = R(r).cfg; cfg.rec.prior = pr{1};
            cfg.tag = ['R2_Y5_' pr{1} '_' R(r).name]; cfg.notes = 'R2 Y5: prior, all baselines';
            out = kt_exp_decim(cfg); kt_write_rows(fid, cfg.tag, R(r).name, out.rows);
        end
    end
    fclose(fid);
end

%% Y6: HFOV with the duty bounded by 1/D (m3)
if any(strcmp(RUN, 'Y6'))
    fid = open_csv('Y6_hfov_duty');
    r = find(strcmp({R.name}, 'hfov'));
    cfg = R(r).cfg; cfg.decim.duty_max = 1;
    cfg.tag = 'R2_Y6_hfov_dutymax1'; cfg.notes = 'R2 Y6: duty <= 1/D at the reduced rate';
    out = kt_exp_decim(cfg); kt_write_rows(fid, cfg.tag, 'hfov', out.rows);
    fclose(fid);
end

%% Y7: stagger inside an injection, native frame rate (m1)
if any(strcmp(RUN, 'Y7'))
    fid = open_csv('Y7_stagger');
    for r = 1:numel(R)
        cfg = R(r).cfg; cfg.stagger.kappa = 0:0.25:1;
        cfg.tag = ['R2_Y7_' R(r).name]; cfg.notes = 'R2 Y7: per-injection vs staggered timing';
        out = kt_stagger_test(cfg); kt_write_rows(fid, cfg.tag, R(r).name, out.rows);
    end
    fclose(fid);
end

%% Y8: simulation with model misspecification (m12), linear forward model, 30 dB
if any(strcmp(RUN, 'Y8'))
    s0 = kt_config();
    s0.sim.forward = 'linear'; s0.sim.snr_db = 30;
    s0.rec.timing_source = 'estimate'; s0.rec.timing_method = 'map';
    s0.rec.timing_assumed.mode = 'per_stim';
    kt_init_eidors(s0);
    Mdl = kt_build_models(s0, true);
    T = 1/s0.acq.fs;
    % name, acq.timing, jitter (s), f0 drift (relative change over the record)
    S = { 'reference',   'per_stim', 0,          0;
          'stagger',     'per_meas', 0,          0;
          'jitter25',    'per_stim', 0.25*T/16,  0;
          'jitter50',    'per_stim', 0.50*T/16,  0;
          'drift0p2',    'per_stim', 0,          0.002;
          'drift0p5',    'per_stim', 0,          0.005;
          'combined',    'per_meas', 0.25*T/16,  0.002 };
    fid = open_csv('Y8_misspec');
    for k = 1:size(S, 1)
        for sd = 1:20
            cfg = s0; cfg.seed = sd;
            cfg.acq.timing = S{k, 2}; cfg.sim.jitter_s = S{k, 3}; cfg.sim.f0_drift = S{k, 4};
            out = kt_sim_core(cfg, Mdl);
            kt_write_rows(fid, sprintf('R2_Y8_%s_seed%d', S{k, 1}, sd), S{k, 1}, out.rows);
        end
        fprintf('Y8 %s done\n', S{k, 1});
    end
    fclose(fid);
end
