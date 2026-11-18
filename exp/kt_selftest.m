function ok = kt_selftest()
%KT_SELFTEST Checks of the k-t core on a synthetic linear problem (no EIDORS).
%   ok = kt_selftest()   run this first; every line must print PASS.
%   Problem: 16 injections x 13 measurements, ~300 pixels in a unit disk,
%   smooth random Jacobian, fs = 44 Hz, components at 10/20/30 Hz (HFOV),
%   2 Hz (cardiac); a second case puts f0 = 11 Hz so that the 3rd harmonic
%   (33 Hz) folds onto the fundamental's bin.
%
% Author: Hoang Nhut Huynh

rng(1);
fs = 44; S = 16; K = 13; F = 440;
stim_idx = repelem((1:S).', K);
k_idx = repmat((1:K).', S, 1);
M = S*K;
[gx, gy] = meshgrid(linspace(-1, 1, 20));
in = gx.^2 + gy.^2 <= 1;
cent = [gx(in), gy(in)];
N = size(cent, 1);
area = ones(N, 1);
Ksm = exp(-pdist2_(cent, cent).^2 / (2*0.25^2));
J = randn(M, N) * Ksm / N;

P = struct('n_stim', S, 'mode', 'per_stim', 'duty', 0.9, 'shift', 3, 'dir', 1);
tau = kt_timing_model(stim_idx, k_idx, fs, P);
lungs = (abs(cent(:,1)) > 0.2) & (cent(:,1).^2 + cent(:,2).^2 < 0.8);
rlung = lungs & cent(:,1) < 0;
heart = (cent(:,1).^2 + (cent(:,2)+0.4).^2) < 0.04;

ok = true;

% --- alias map
C = kt_alias_map([10 20 30 40], fs, 0.1);
ok = check('alias map bins/signed/ell', isequal([C.bin], [10 20 14 4]) && ...
    isequal([C.signed], [10 20 -30 -40]) && isequal([C.ell], [0 0 -1 -1]), ok);

% --- case 1: HFOV 10 Hz + harmonics, no collision
nu = [10 20 30 2];
X = [-0.05*lungs, -0.012*exp(1i*0.7)*lungs, -0.01*exp(1i*1.4)*rlung, 0.01*heart];
[dv_seq, dv_sim] = acquire(J, X, nu, tau, fs, F, 80);
w = kt_noise_floor(dv_seq, fs, [10 20 14 2], 'hann', 4*fs/F);
rec = kt_recon_matrix(J, w, 'noser', 0.03);
opt = struct('window', 'hann', 'tol', 2*fs/F, 'l_max', 2);
meth = {'naive', 'yerworth', 'kt_single', 'kt', 'kt_auto'};
A = kt_analyze(dv_seq, fs, nu, tau, rec, meth, opt);
R0 = kt_analyze(dv_sim, fs, nu, tau, rec, {'naive'}, opt);
rho = zeros(numel(meth), numel(nu));
for i = 1:numel(meth)
    for c = 1:numel(nu)
        m = kt_compare(A.X.(meth{i})(:, c), R0.X.naive(:, c), area);
        rho(i, c) = m.rho;
    end
end
fprintf('rho vs simultaneous-acquisition reference (rows %s; cols 10,20,30,2 Hz)\n', strjoin(meth, ', '));
disp(round(rho*1e4)/1e4);
ok = check('kt equals simultaneous reference (rho > 0.999)', all(rho(4, :) > 0.999), ok);
ok = check('yerworth fails on the aliased 30 Hz harmonic (rho < 0.95)', rho(2, 3) < 0.95, ok);
ok = check('kt_auto picks the physical alias of every component', ...
    isequal(A.nu_used.kt_auto, [C(1:3).signed, 2]), ok);

% --- timing estimation from the 10 Hz data
y0 = kt_spectrum(dv_seq, fs, 10, 'hann');
[Pe, ti] = kt_estimate_timing(y0, 10, stim_idx, k_idx, fs, P, 0:0.01:1, w);
fprintf('timing estimate: duty %.2f shift %d dir %d (true 0.90 3 1), coherence %.4f (no timing %.4f)\n', ...
    Pe.duty, Pe.shift, Pe.dir, ti.c_best, ti.c_no_timing);
ok = check('timing estimate recovers duty/shift/dir', abs(Pe.duty - 0.9) <= 0.02 && ...
    Pe.shift == 3 && Pe.dir == 1, ok);

% --- case 2: f0 = 11 Hz, 33 Hz folds onto the 11 Hz bin
nu2 = [11 33 2];
X2 = [-0.05*lungs, -0.02*exp(1i*0.9)*rlung, 0.01*heart];
[dv2, ~] = acquire(J, X2, nu2, tau, fs, F, 60);
w2 = kt_noise_floor(dv2, fs, [11 2], 'hann', 4*fs/F);
rec2 = kt_recon_matrix(J, w2, 'noser', 0.03);
C2 = kt_alias_map(nu2, fs, opt.tol);
ok = check('collision detected (11 and 33 Hz in one group)', C2(1).group == C2(2).group, ok);
A2 = kt_analyze(dv2, fs, nu2, tau, rec2, {'kt_single', 'kt'}, opt);
r_single = kt_compare(A2.X.kt_single(:, 2), X2(:, 2), area);
r_joint = kt_compare(A2.X.kt(:, 2), X2(:, 2), area);
fprintf('33 Hz vs truth: kt_single rho %.3f amp %.2f | kt joint rho %.3f amp %.2f\n', ...
    r_single.rho, r_single.amp, r_joint.rho, r_joint.amp);
ok = check('joint solve separates the colliding harmonic better', r_joint.rho > r_single.rho, ok);
sep = kt_separability(J, w2, tau, -44, [10 50 100]);
sep0 = kt_separability(J, w2, zeros(size(tau)), -44, [10 50 100]);
fprintf('min principal angle (deg) r=10/50/100: sequential %s | simultaneous %s\n', ...
    mat2str(round(sep.min_angle_deg*10)/10), mat2str(round(sep0.min_angle_deg*10)/10));
ok = check('simultaneous sampling cannot separate the aliases (angle 0)', all(sep0.min_angle_deg < 1e-3), ok);

if ok
    fprintf('\nkt_selftest: ALL PASS\n');
else
    fprintf('\nkt_selftest: SOME CHECKS FAILED\n');
end
end

function [dv_seq, dv_sim] = acquire(J, X, nu, tau, fs, F, snr_db)
tf = (0:F-1)/fs;
P = J*X;
dv_seq = zeros(size(J, 1), F);
dv_sim = dv_seq;
for c = 1:numel(nu)
    dv_seq = dv_seq + real(P(:, c) .* exp(1i*2*pi*nu(c)*(tau(:) + tf)));
    dv_sim = dv_sim + real(P(:, c) .* exp(1i*2*pi*nu(c)*(0*tau(:) + tf)));
end
sd = sqrt(mean(dv_seq(:).^2)) / 10^(snr_db/20);
dv_seq = dv_seq + sd*randn(size(dv_seq));
dv_sim = dv_sim + sd*randn(size(dv_sim));
end

function ok = check(name, cond, ok)
if cond
    fprintf('PASS  %s\n', name);
else
    fprintf('FAIL  %s\n', name);
    ok = false;
end
end

function D = pdist2_(a, b)
D = sqrt(max(sum(a.^2, 2) + sum(b.^2, 2).' - 2*a*b.', 0));
end
