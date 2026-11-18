% run_R1_revision  Analyses requested in the simulated review (round R1) of the k-t EIT manuscript.
% Sections can be run separately (Ctrl+Enter on a %% cell) or all at once (about 1.5 to 2 h,
% most of it in X2). Metrics only, no figures. Each section writes to runs/ as usual; X2 also
% writes one table runs/R1_X2_nonlinear_sim.csv.
%
%   X2  nonlinear forward model (EIDORS fwd_solve on the fine mesh), 20 noise seeds x SNR 20/30/40 dB,
%       for the duty-mismatch sweep (Fig. 2b) and the f0 sweep (Fig. 2c)
%   X5  sensitivity of the reduced-rate reconstruction to an error in f0 (surfactant record)
%   X6  physiological bound on the alias order chosen by k-t (auto) for a cardiac line
%   X7  EIDORS frame_timing_correct in 'linear' (interpolation) mode as a second FTC baseline
%   X4  two further reconstruction matrices (Laplace and Tikhonov priors) for generality
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;
if ~exist('RUN', 'var'), RUN = {'X2', 'X5', 'X6', 'X7', 'X4'}; end   % e.g. RUN = {'X7'}; run_R1_revision

%% X2: nonlinear forward model, noise ensemble
if any(strcmp(RUN, 'X2'))
    base = kt_config();
    base.sim.forward = 'nonlinear';
    kt_init_eidors(base);
    Mdl = kt_build_models(base, true);
    seeds = 1:20; snrs = [20 30 40];
    sweeps = {'rec.timing_assumed.duty', [0.5 0.7 0.8 0.85 0.9 0.95 1.0];
              'sim.f0', 6:1:20};
    fid = fopen(fullfile('runs', 'R1_X2_nonlinear_sim.csv'), 'w');
    fprintf(fid, 'sweep,value,snr_db,seed,component,method,metric,metric_value\n');
    t0 = tic;
    for s = 1:size(sweeps, 1)
        for v = sweeps{s, 2}
            for snr = snrs
                for sd = seeds
                    cfg = base; cfg.seed = sd; cfg.sim.snr_db = snr;
                    if s == 1, cfg.rec.timing_source = 'assumed'; end
                    cfg = kt_set(cfg, sweeps{s, 1}, v);
                    out = kt_sim_core(cfg, Mdl);
                    r = out.rows;
                    for i = 1:size(r, 1)
                        if any(strcmp(r{i, 5}, {'rho_true', 'relerr_true', 'rho_ref'}))
                            fprintf(fid, '%s,%g,%g,%d,%s,%s,%s,%.6g\n', sweeps{s, 1}, v, snr, sd, r{i, 3}, r{i, 4}, r{i, 5}, r{i, 6});
                        end
                    end
                end
                fprintf('X2 %s = %g, SNR %g dB done (%.0f s)\n', sweeps{s, 1}, v, snr, toc(t0));
            end
        end
    end
    fclose(fid);
end

%% common settings for the frame-removal analyses (X4 to X7)
dbase = kt_config();
dbase.save_figs = false; dbase.save_mat = false;
dbase.rec.eidors_temporal = {};
dbase.model.shape = 'pig_23kg';
surf = dbase;                                  % PEEP-trial pig after surfactant
surf.real.files = {fullfile(pwd, 'Data', 'pig_2003', 'p1130122.get')};
surf.real.fs = 13; surf.real.f0_search = [0.6 0.8]; surf.real.n_harm = 9;
piglet = dbase;                                % piglet 2, cardiac line
piglet.real.files = {fullfile(pwd, 'Data', 'Pig_torso', 'g10-c005.get')};
piglet.real.fs = 13; piglet.real.f0_search = [0.30 0.42]; piglet.real.n_harm = 9;
piglet.real.cardiac_search = [3.40 3.65];

%% X5: error in f0 (only the reduced-rate reconstruction is perturbed; the reference is not)
if any(strcmp(RUN, 'X5'))
    for e = [-0.01 -0.005 -0.002 -0.001 0 0.001 0.002 0.005 0.01]
        cfg = surf; cfg.decim.f0_err = e;
        cfg.rec.methods = {'naive', 'yerworth', 'kt'}; cfg.rec.try_eidors_ftc = false;
        cfg.tag = sprintf('R1_X5_surf_f0err_%+.1fpct', 100*e);
        cfg.notes = 'R1 X5: f0 error at the reduced rate';
        kt_exp_decim(cfg);
    end
end

%% X6: physiological bound for k-t (auto) on the cardiac line (|nu| <= 5 Hz)
if any(strcmp(RUN, 'X6'))
    cfg = kt_config();
    cfg.real.files = {fullfile(pwd, 'Data', 'chb_hfov_porcine', 'CHB_pig06-55.get')};
    cfg.real.cardiac_search = [1.5 3.5];
    cfg.model.shape = 'pig_23kg';
    cfg.rec.timing_source = 'estimate'; cfg.rec.timing_method = 'map'; cfg.real.cv_timing = false;
    cfg.rec.try_eidors_ftc = false;
    for nm = [inf 5]
        cfg.rec.cardiac_numax = nm;
        cfg.tag = sprintf('R1_X6_hfov_numax%g', nm);
        cfg.notes = 'R1 X6: kt_auto with a physiological bound on the cardiac line';
        kt_exp_real(cfg);
    end
    cfg = kt_config(); cfg.sim.f0 = 14;            % simulation: 3rd harmonic (42 Hz) on the 2 Hz cardiac line
    for nm = [inf 5]
        cfg.rec.cardiac_numax = nm;
        cfg.tag = sprintf('R1_X6_sim_f0_14_numax%g', nm);
        cfg.notes = 'R1 X6: kt_auto bound, simulation f0 = 14 Hz';
        kt_exp_sim(cfg);
    end
end

%% X7: FTC in interpolation mode
if any(strcmp(RUN, 'X7'))
    for R = {surf, piglet; 'surf', 'piglet2'}
        cfg = R{1}; cfg.rec.ftc_method = 'linear';
        cfg.rec.methods = {'naive', 'yerworth', 'kt'};
        cfg.tag = ['R1_X7_' R{2} '_ftc_linear'];
        cfg.notes = 'R1 X7: frame_timing_correct linear mode';
        kt_exp_decim(cfg);
    end
end

%% X4: other reconstruction matrices (priors)
if any(strcmp(RUN, 'X4'))
    for R = {surf, piglet; 'surf', 'piglet2'}
        for pr = {'laplace', 'tikhonov'}
            cfg = R{1}; cfg.rec.prior = pr{1};
            cfg.rec.methods = {'naive', 'yerworth', 'kt'};
            cfg.tag = ['R1_X4_' R{2} '_' pr{1}];
            cfg.notes = 'R1 X4: reconstruction matrix with another prior';
            kt_exp_decim(cfg);
        end
    end
end
