function [P, info] = kt_estimate_timing(y0, f0, stim_idx, k_idx, fs, P0, duty_grid, w)
%KT_ESTIMATE_TIMING Acquisition schedule estimated from the data at one frequency.
%   [P, info] = kt_estimate_timing(y0, f0, stim_idx, k_idx, fs, P0, duty_grid, w)
%   y0 : M x 1 data at an unaliased frequency f0 (< fs/2) with a strong,
%        spatially in-phase oscillation (e.g. the HFOV fundamental).
%   If the oscillation is in phase, J*X has one common phase, so
%   z_m = y0_m * exp(-1i*2*pi*f0*tau_m) is real up to a global phase and the
%   coherence  c = |sum w z.^2| / sum w |z|.^2  reaches 1 at the true schedule.
%   Grid search over dir (+1,-1), shift (0..S-1) and duty (duty_grid >= 0).
%   P0 supplies n_stim and mode. info.c(dir, shift+1, duty) is the surface,
%   info.c_no_timing the coherence with timing ignored.
%   Caveat: phase-heterogeneous oscillation (pendelluft) lowers c_best; the
%   arg max is then less reliable. Report c_best in every run.
%
% Author: Hoang Nhut Huynh

if nargin < 8 || isempty(w), w = ones(size(y0)); end
y0 = y0(:);
w = w(:);
S = P0.n_stim;
dirs = [1 -1];
cs = zeros(2, S, numel(duty_grid));
P = P0;
for a = 1:2
    for s = 0:S-1
        for d = 1:numel(duty_grid)
            P.dir = dirs(a); P.shift = s; P.duty = duty_grid(d);
            cs(a, s+1, d) = coherence(y0, w, f0, kt_timing_model(stim_idx, k_idx, fs, P));
        end
    end
end
[cmax, imax] = max(cs(:));
[a, s, d] = ind2sub(size(cs), imax);
P.dir = dirs(a); P.shift = s - 1; P.duty = duty_grid(d);
info.c = cs;
info.c_best = cmax;
info.duty_grid = duty_grid;
info.c_no_timing = coherence(y0, w, f0, zeros(size(y0)));
end

function c = coherence(y0, w, f0, tau)
z = y0 .* exp(-1i*2*pi*f0*tau);
c = abs(sum(w .* z.^2)) / sum(w .* abs(z).^2);
end
