function rows = kt_order_design(cfg, Mdl, order)
%KT_ORDER_DESIGN Fixed versus per-frame random injection order (acquisition design, R3 review).
%   rows = kt_order_design(cfg, Mdl, order)
%   order = 'fixed'  : equal-slot schedule of cfg.acq ('per_stim'), the same in every frame
%         = 'random' : the S injection slots are permuted independently in every frame;
%                      the schedule is taken as known (logged by the device)
%   Linear or nonlinear forward model as in cfg.sim; the schedule is known in both cases.
%   Methods ('<order>_<method>' in the output):
%     fixed  : 'yerworth', 'kt_single', 'kt' (joint)            via kt_analyze
%              'harmreg'  per-channel least-squares fit of all component frequencies on the
%                         acquisition times t_{f,m}, then the shared one-step reconstructor R
%                         (for an isolated component this equals kt_single; for colliding
%                         components the per-channel fit is rank deficient)
%     random : 'nudft'    per-channel Hann-windowed DFT at the true frequency on t_{f,m}, then R
%              'harmreg'  as above (all components fitted jointly in every channel)
%   Scores: rho_true, relerr_true (phantom on the inverse mesh) and rho_ref (reconstruction
%   of the same phantom acquired simultaneously). rows: {param, pvalue, component, method,
%   metric, value} with param = 'f0'.
%
% Author: Hoang Nhut Huynh

assert(strcmp(cfg.acq.timing, 'per_stim'), 'kt_order_design: cfg.acq.timing must be ''per_stim''');
fs = cfg.acq.fs; T = 1/fs;
I = Mdl.inv;
F = round(cfg.sim.duration * fs);
tf = (0:F-1) / fs;
S = I.n_stim;
P = struct('n_stim', S, 'mode', 'per_stim', 'duty', cfg.acq.duty, 'shift', cfg.acq.shift, 'dir', cfg.acq.dir);
switch order
    case 'fixed'
        Tf = kt_timing_model(Mdl.fwd.stim_idx, Mdl.fwd.k_idx, fs, P);
        Ti = kt_timing_model(I.stim_idx, I.k_idx, fs, P);
    case 'random'
        s_old = rng; rng(cfg.seed + 104729);
        slot = zeros(S, F);
        for f = 1:F, slot(:, f) = randperm(S).' - 1; end
        rng(s_old);
        Tf = (cfg.acq.duty * T / S) * slot(Mdl.fwd.stim_idx(:), :);
        Ti = (cfg.acq.duty * T / S) * slot(I.stim_idx(:), :);
    otherwise
        error('kt_order_design: unknown order "%s"', order);
end
[dv_seq, dv_sim] = kt_simulate(cfg, Mdl, Tf);

[nu_all, Xt_all, names_all, analyze] = kt_phantom(cfg, I.cent, I.R);
ia = find(analyze);
nu = nu_all(ia); Xt = Xt_all(:, ia); names = names_all(ia);
opt = struct('window', cfg.rec.window, 'tol', 2*fs/F, 'l_max', cfg.rec.l_max);
% Noise weights from the per-channel residual of the harmonic fit (all components modeled).
% v1 used kt_noise_floor, whose signal-free frame bins contain the aliasing leakage of a
% random order, so the weights followed the signal, not the noise, and degraded with SNR.
% Only relative weights matter for the NOSER prior (lam2 = alpha), so time-domain variances suffice.
Tm = Ti + tf;                                   % M x F acquisition times (inverse rows)
[A, s2] = harm_regress(dv_seq, Tm, nu_all);
w = 1 ./ s2;
prior = cfg.rec.prior;
if strcmp(prior, 'laplace'), prior = 'external'; end
rec = kt_recon_matrix(I.J, w, prior, cfg.rec.alpha, Mdl.RtR_laplace);

R0 = kt_analyze(dv_sim, fs, nu, zeros(size(I.J, 1), 1), rec, {'naive'}, opt);
Xsim = R0.X.naive;

X = struct();
valid = true(1, numel(nu));
X.harmreg = rec.R * A(:, ia);
switch order
    case 'fixed'
        Ak = kt_analyze(dv_seq, fs, nu, Ti, rec, {'yerworth', 'kt_single', 'kt'}, opt);
        X.yerworth = Ak.X.yerworth; X.kt_single = Ak.X.kt_single; X.kt = Ak.X.kt;
        valid = [Ak.C.valid];
    case 'random'
        X.nudft = rec.R * nudft(dv_seq, Tm, nu, cfg.rec.window);
end

rows = cell(0, 6);
mlist = fieldnames(X);
for c = 1:numel(nu)
    if ~valid(c), continue; end
    for k = 1:numel(mlist)
        x = X.(mlist{k})(:, c);
        st = kt_compare(x, Xt(:, c), I.area);
        sr = kt_compare(x, Xsim(:, c), I.area);
        mname = [order '_' mlist{k}];
        rows = [rows; {'f0', cfg.sim.f0, names{c}, mname, 'rho_true', st.rho; ...
                       'f0', cfg.sim.f0, names{c}, mname, 'relerr_true', st.relerr; ...
                       'f0', cfg.sim.f0, names{c}, mname, 'rho_ref', sr.rho}]; %#ok<AGROW>
    end
    if strcmp(order, 'fixed')         % simultaneous reference itself, for scale
        st = kt_compare(Xsim(:, c), Xt(:, c), I.area);
        rows = [rows; {'f0', cfg.sim.f0, names{c}, 'simref', 'rho_true', st.rho}]; %#ok<AGROW>
    end
end
end

function [A, s2] = harm_regress(dv, t, nu)
% Per-channel least squares: dv_m(t) = a0 + sum_k Re{a_k exp(i 2 pi nu_k t)}; returns M x K a_k
% and the residual variance s2 (M x 1) of every channel.
[M, F] = size(dv);
K = numel(nu);
A = zeros(M, K);
s2 = zeros(M, 1);
for m = 1:M
    tm = t(m, :).';
    B = [ones(F, 1), cos(2*pi*tm*nu(:).'), sin(2*pi*tm*nu(:).')];
    c = lsqminnorm(B, dv(m, :).');
    A(m, :) = (c(2:K+1) - 1i*c(K+2:end)).';
    s2(m) = sum((dv(m, :).' - B*c).^2) / (F - numel(c));
end
end

function Y = nudft(dv, t, nu, win)
% Windowed DFT of every channel at the true frequencies on its own acquisition times
% (same scaling as kt_spectrum: a*cos(2*pi*nu*t + p) gives a*exp(1i*p)).
F = size(dv, 2);
n = (0:F-1);
switch win
    case 'hann', wv = 0.5 - 0.5*cos(2*pi*(n + 0.5)/F);
    otherwise,   wv = ones(1, F);
end
dv = dv - mean(dv, 2);
Y = zeros(size(dv, 1), numel(nu));
for k = 1:numel(nu)
    Y(:, k) = 2 * sum(dv .* wv .* exp(-1i*2*pi*nu(k)*t), 2) / sum(wv);
end
end
