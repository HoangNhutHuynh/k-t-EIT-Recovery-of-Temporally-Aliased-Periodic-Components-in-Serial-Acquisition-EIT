% run_R3_revision  Analyses for the full review (R3: EIC + 3 reviewers + Devil's Advocate).
% Metrics only, no figures. Every section appends its rows to one CSV in runs/:
%   runs/R3_<section>.csv  columns: tag, record, param, pvalue, component, method, metric, value
% Approximate run times (desktop PC): Z1 35 min, Z2 15 min.
%
%   Z1 (risk 4)  blind model choice for the joint solve: for every aliased component the
%                k-t model with or without an unaliased line at its folded bin is chosen by
%                the negative log evidence (method 'kt_sel'), nine records, alpha = 0.03
%   Z2 (risk 1, R3-1, R1-1)  acquisition design in simulation: fixed injection order (k-t,
%                joint k-t, per-channel harmonic regression) against an injection order
%                permuted in every frame (per-channel DFT and harmonic regression on the
%                known acquisition times); f0 = 10 (isolated aliased h3), 11 (h3 on the
%                fundamental) and 14 Hz (h3 on the cardiac line); 20, 30, 40 dB; 20 seeds
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;
if ~exist('RUN', 'var'), RUN = {'Z1', 'Z2'}; end   % e.g. RUN = {'Z2'}; run_R3_revision
if ~exist('runs', 'dir'), mkdir('runs'); end   % fresh clone: runs/ is not tracked
open_csv = @(s) fopen(fullfile('runs', ['R3_' s '.csv']), 'a');

%% Z1: blind (evidence-based) choice of the nuisance line, nine records
if any(strcmp(RUN, 'Z1'))
    base = kt_config();
    base.save_figs = false; base.save_mat = false;
    base.rec.eidors_temporal = {};
    base.real.extra_search = {};
    R = kt_r2_records(base);
    fid = open_csv('Z1_evidence');
    for r = 1:numel(R)
        cfg = R(r).cfg; cfg.decim.nuisance_select = 'evidence';
        cfg.rec.methods = {'yerworth', 'kt_single', 'kt'}; cfg.rec.try_eidors_ftc = false; cfg.rec.eidors_temporal = {};
        cfg.tag = ['R3_Z1_' R(r).name]; cfg.notes = 'R3 Z1: evidence-based choice of the nuisance line';
        out = kt_exp_decim(cfg); kt_write_rows(fid, cfg.tag, R(r).name, out.rows);
    end
    fclose(fid);
end

%% Z2: fixed versus per-frame random injection order (simulation, linear forward model)
if any(strcmp(RUN, 'Z2'))
    s0 = kt_config();
    s0.save_figs = false; s0.save_mat = false;
    s0.sim.forward = 'linear'; s0.acq.timing = 'per_stim';
    kt_init_eidors(s0);
    Mdl = kt_build_models(s0, true);
    fid = open_csv('Z2_order');
    for f0 = [10 11 14]
        c0 = s0; c0.sim.f0 = f0;
        if f0 == 11                % as E3: h2 at Nyquist removed, h3 on the fundamental
            keep = ismember({c0.sim.components.name}, {'hfov_f0', 'hfov_h3', 'cardiac', 'drift'});
            c0.sim.components = c0.sim.components(keep);
            c0.sim.components(strcmp({c0.sim.components.name}, 'hfov_h3')).amp = -0.02;
        end
        for snr = [20 30 40]
            for sd = 1:20
                cfg = c0; cfg.seed = sd; cfg.sim.snr_db = snr;
                for order = {'fixed', 'random'}
                    rows = kt_order_design(cfg, Mdl, order{1});
                    kt_write_rows(fid, sprintf('R3_Z2_f%d_snr%d_seed%d', f0, snr, sd), sprintf('snr%d', snr), rows);
                end
            end
            fprintf('Z2 f0 %d Hz, %d dB done\n', f0, snr);
        end
    end
    fclose(fid);
end
