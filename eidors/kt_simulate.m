function [dv_seq, dv_sim] = kt_simulate(cfg, Mdl, tau)
%KT_SIMULATE Sequential and simultaneous acquisition of the phantom time series.
%   [dv_seq, dv_sim] = kt_simulate(cfg, Mdl, tau)
%   dv = v./vh - 1 on the fine forward mesh (M x F), with white noise at
%   cfg.sim.snr_db. dv_seq: measurement m of frame f taken at f/fs + tau(m);
%   dv_sim: all measurements at f/fs (ideal reference, not realisable).
%   cfg.sim.forward = 'linear'    Jacobian of the fine mesh (fast)
%                   = 'nonlinear' one fwd_solve per distinct acquisition time
%                                 (16 per frame for 'per_stim'; slow for 'per_meas')
%
% Author: Hoang Nhut Huynh

persistent cache
fs = cfg.acq.fs;
F = round(cfg.sim.duration * fs);
tf = (0:F-1) / fs;
I = Mdl.fwd;
[nu, X] = kt_phantom(cfg, I.cent, I.R);
if size(tau, 2) > 1        % R3: per-frame schedule (M x F), e.g. injection order permuted per frame
    assert(size(tau, 2) == F, 'kt_simulate: tau must be M x 1 or M x F');
    t_seq = tau + tf;
else
    t_seq = tau(:) + tf;
end
t_sim = zeros(size(tau, 1), 1) + tf;
% model misspecification (R2 review, linear forward only):
%   sim.jitter_s : sd (s) of a random offset of every measurement in every frame
%   sim.f0_drift : relative linear change of the harmonic frequencies over the record
%                  (instantaneous frequency nu*(1 + drift*(t/T_rec - 1/2)))
jit = kt_get(cfg, 'sim.jitter_s', 0);
drift = kt_get(cfg, 'sim.f0_drift', 0);
if (jit > 0 || drift ~= 0) && ~strcmp(cfg.sim.forward, 'linear')
    error('kt_simulate: sim.jitter_s and sim.f0_drift need sim.forward = ''linear''');
end
if jit > 0
    s_old = rng; rng(cfg.seed + 7919);
    t_seq = t_seq + jit*randn(size(t_seq));
    rng(s_old);
end
harm = [cfg.sim.components.harm] > 0;
Trec = cfg.sim.duration;
warp = @(t, c) t + (drift*harm(c)) * (t.^2/(2*Trec) - t/2);   % phase time of component c
switch cfg.sim.forward
    case 'linear'
        P = I.J * X;
        dv_seq = zeros(size(tau, 1), F);
        dv_sim = dv_seq;
        for c = 1:numel(nu)
            dv_seq = dv_seq + real(P(:, c) .* exp(1i*2*pi*nu(c)*warp(t_seq, c)));
            dv_sim = dv_sim + real(P(:, c) .* exp(1i*2*pi*nu(c)*warp(t_sim, c)));
        end
    case 'nonlinear'
        % the noise-free data depend on the phantom and the schedule only, so they are
        % cached across calls that change just the noise (seed, SNR) or the reconstruction
        wv = (1:size(X, 1)).';
        key = sprintf('%d|%d|%.6g|%s|%.12g|%.12g|%.12g|%.12g', numel(tau), F, fs, sprintf('%.9g,', nu), ...
                      sum(real(X) .* wv, 'all'), sum(imag(X) .* wv, 'all'), sum(tau(:) .* (1:numel(tau)).'), sum(tau(:)));
        if isempty(cache), cache = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
        if isKey(cache, key)
            c = cache(key); dv_seq = c{1}; dv_sim = c{2};
        else
            dv_seq = acquire_nl(I, X, nu, t_seq);
            dv_sim = acquire_nl(I, X, nu, t_sim);
            if cache.Count > 40, remove(cache, keys(cache)); end
            cache(key) = {dv_seq, dv_sim};
        end
    otherwise
        error('kt_simulate: unknown forward mode "%s"', cfg.sim.forward);
end
rng(cfg.seed);
sd = sqrt(mean(dv_seq(:).^2)) / 10^(cfg.sim.snr_db/20);
dv_seq = dv_seq + sd*randn(size(dv_seq));
dv_sim = dv_sim + sd*randn(size(dv_sim));
end

function dv = acquire_nl(I, X, nu, t)
[M, F] = size(t);
dv = zeros(M, F);
img = mk_image(I.fmdl, 1);
tic_last = tic;
for f = 1:F
    ts = unique(t(:, f));
    for k = 1:numel(ts)
        rows = t(:, f) == ts(k);
        img.elem_data = 1 + real(X * exp(1i*2*pi*nu(:)*ts(k)));
        v = fwd_solve(img);
        dv(rows, f) = v.meas(rows) ./ I.vh(rows) - 1;
    end
    if toc(tic_last) > 20
        fprintf('  nonlinear forward: frame %d / %d\n', f, F);
        tic_last = tic;
    end
end
end
