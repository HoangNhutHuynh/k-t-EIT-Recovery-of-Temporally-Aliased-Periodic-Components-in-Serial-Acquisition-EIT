function cfg = kt_config()
%KT_CONFIG Default configuration of the k-t EIT experiments (Idea 1).
%   cfg = kt_config();  then override fields in the run_E*.m scripts, e.g.
%   cfg.tag = 'alpha0p1'; cfg.rec.alpha = 0.1; cfg.notes = 'why this run';
%   Every run stores the full cfg (config.txt) in its run folder, so change
%   values in the scripts rather than here when you try variants.
%
% Author: Hoang Nhut Huynh

cfg.tag   = 'default';      % short label of the trial (file names, logs)
cfg.notes = '';             % hypothesis / reason for this run (goes to the logs)
cfg.seed  = 1;
cfg.root  = fileparts(mfilename('fullpath'));
cfg.out_root = fullfile(cfg.root, 'runs');
cfg.eidors_startup = '';    % e.g. '<path to>/eidors-v3.12-ng/eidors/startup.m' (empty: EIDORS already on the path)
cfg.keep_figs_open = false;
cfg.save_figs = false;       % false while tuning: no figures (metrics.csv is the record); true for final runs
cfg.save_mat  = false;       % results.mat (images, ~4 MB per run) only when needed
cfg.log_keys = {'acq.duty', 'acq.timing', 'rec.timing_source', 'rec.prior', ...
                'rec.alpha', 'sim.forward', 'sim.snr_db', 'sim.f0'};

% ---- acquisition schedule (simulation truth; real data: initial guess)
cfg.acq.fs     = 44;          % frame rate (Hz)
cfg.acq.n_elec = 16;
cfg.acq.timing = 'per_stim';  % 'per_stim' | 'per_meas'
cfg.acq.duty   = 0.9;         % fraction of the frame period used by the injections
cfg.acq.shift  = 0;           % injection slot read first
cfg.acq.dir    = 1;           % +1 / -1 order of the slots

% ---- reconstruction (identical for every method)
cfg.rec.prior  = 'noser';     % 'noser' | 'tikhonov' | 'laplace'
cfg.rec.alpha  = 0.03;        % lam2 = alpha*trace(J'WJ)/trace(RtR)
cfg.rec.window = 'hann';
cfg.rec.methods = {'naive', 'yerworth', 'kt_single', 'kt', 'kt_auto'};
cfg.rec.l_max  = 2;           % alias orders searched by kt_auto
cfg.rec.tol_factor = 2;       % bins closer than tol_factor*fs/F are one bin (joint solve); larger = also near-collisions
cfg.rec.timing_source = 'true';   % 'true' (sim oracle) | 'assumed' | 'estimate'
cfg.rec.timing_assumed = struct('mode', 'per_stim', 'duty', 0.9, 'shift', 0, 'dir', 1);
cfg.rec.timing_ref = 'hfov_f0';   % component used by 'estimate' (must be unaliased)
cfg.rec.timing_method = 'map';    % 'map' (min MAP objective, equal slots) | 'free' (real data: one offset per injection, refined from 'map') | 'coherence' (in-phase f0 only)
cfg.rec.duty_grid = 0:0.01:1;
cfg.rec.try_eidors_ftc = true;    % also try EIDORS frame_timing_correct as a baseline
cfg.rec.eidors_temporal = {};
cfg.rec.cardiac_numax = inf;   % kt_auto: largest |nu| (Hz) allowed for a component named 'cardiac' (physiological bound; inf = none)
cfg.rec.ftc_method = 'FFT';     % EIDORS frame_timing_correct mode: 'FFT' | 'linear'     % real data: extra baselines: 'time_prior' (own, Adler 2007 data form, same J/prior/lambda), 'kalman' (EIDORS, not tuned)

% ---- EIDORS models (2-D circular, built-in, no Netgen needed)
cfg.model.fwd = 'f2c2';       % fine mesh used to simulate data
cfg.model.inv = 'c2c2';       % coarse mesh used to reconstruct
cfg.model.shape = 'circle';   % real data: shape_library outline for the inverse mesh, e.g. 'pig_23kg', 'neonate'
cfg.model.rotate_deg = 0;     % extra rotation of the outline (degrees), check fig_model.png

% ---- simulation
cfg.sim.duration = 10;        % s
cfg.sim.snr_db   = 40;        % w.r.t. rms of the difference signal
cfg.sim.forward  = 'linear';  % 'linear' (fine-mesh Jacobian, fast) | 'nonlinear' (fwd_solve per injection)
cfg.sim.f0       = 10;        % HFOV fundamental (Hz)
% harm > 0: nu = harm*f0 ; harm = 0: nu as given. analyze = component is imaged and scored.
% phase_lr_deg: extra phase of the image-right half (x > 0, patient left lung), e.g. pendelluft-like lag.
cfg.sim.components = struct( ...
    'name',      {'hfov_f0', 'hfov_h2', 'hfov_h3', 'cardiac', 'drift'}, ...
    'harm',      {1,         2,         3,         0,         0}, ...
    'nu',        {NaN,       NaN,       NaN,       2.0,       0.15}, ...
    'amp',       {-0.05,     -0.012,    -0.010,    0.010,     -0.02}, ...
    'phase_deg', {0,         40,        80,        0,         0}, ...
    'region',    {'lungs',   'lungs',   'right_lung', 'heart', 'lungs'}, ...
    'phase_lr_deg', {0,      0,         0,         0,         0}, ...
    'analyze',   {true,      true,      true,      true,      false});

% ---- real data (chb_hfov_porcine by default)
cfg.real.files     = {fullfile(fileparts(mfilename('fullpath')), 'Data', 'chb_hfov_porcine', 'CHB_pig06-55.get')};  % see README, Data
cfg.real.reader    = 'goe_raw';    % 'goe_raw' (own reader, device order, drops unstored columns) | 'eidors'
cfg.real.formats   = {'GET-RAW', 'GET'};   % tried in this order by eidors_readdata
cfg.real.fs        = 44;        % data page: 44 Hz (.prl header is wrong)
cfg.real.t_range   = [0 Inf];   % s, segment analysed (look at fig_global first)
cfg.real.f0_search = [8 12];    % Hz, window for the HFOV fundamental
cfg.real.n_harm    = 4;
cfg.real.extra_search = {};     % cell of Hz bands; the global-spectrum peak in each is added as an unaliased line (line1, line2, ...)
cfg.real.cardiac_search = [];   % Hz band, e.g. [1.5 3.5]: add the cardiac line found in the global spectrum
cfg.real.extra_fold_of = [];   % harmonics h whose folded bin also gets an unaliased component (joint solve), e.g. 4
cfg.real.extra_excl = [];       % Hz, extra lines kept out of the noise floor (e.g. cardiac)
cfg.real.cv_timing = false;     % two-fold held-out check of the estimated schedule
cfg.real.block_sec = 5;         % s, block length for the block bootstrap
cfg.real.n_boot    = 2000;

% ---- frame-rate reduction test on real data (E7)
cfg.decim.factor = 2;            % keep every factor-th frame
cfg.decim.timing = 'decimated';  % 'decimated' (schedule re-estimated at the low rate) | 'full'
cfg.decim.fig_rows = {};
cfg.decim.f0_err = 0;            % relative error put on f0 for the reduced-rate reconstruction only (sensitivity test, e.g. 0.002)         % component names for an extra main-text figure, e.g. {'h5','h6','h7','h8'}
end
