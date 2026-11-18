function [P, info] = kt_estimate_timing_map(Y, nus, stim_idx, k_idx, fs, P0, duty_grid, rec)
%KT_ESTIMATE_TIMING_MAP Acquisition schedule by minimum MAP objective (no in-phase assumption).
%   [P, info] = kt_estimate_timing_map(Y, nus, stim_idx, k_idx, fs, P0, duty_grid, rec)
%   Y   : M x K data at the bins of K components (kt_spectrum)
%   nus : 1 x K signed physical frequencies of those components (kt_alias_map .signed)
%   For one frequency the normal matrix J'WJ + lam2*RtR does not depend on the
%   schedule (D is unitary and diagonal), so the minimum objective Phi(tau) is
%   -2 log evidence up to a constant: comparable across schedules. The images
%   are complex, so spatially varying phase (e.g. pendelluft) is allowed.
%   Score = sum_k Phi_k; grid over dir (+1,-1), shift (0..S-1), duty_grid.
%   info.phi(dir, shift+1, duty), info.phi_best, info.phi_no_timing (D = I).
%
% Author: Hoang Nhut Huynh

S = P0.n_stim;
dirs = [1 -1];
phi = zeros(2, S, numel(duty_grid));
P = P0;
for a = 1:2
    for s = 0:S-1
        for d = 1:numel(duty_grid)
            P.dir = dirs(a); P.shift = s; P.duty = duty_grid(d);
            phi(a, s+1, d) = score(rec, Y, nus, kt_timing_model(stim_idx, k_idx, fs, P));
        end
    end
end
[pmin, imin] = min(phi(:));
[a, s, d] = ind2sub(size(phi), imin);
P.dir = dirs(a); P.shift = s - 1; P.duty = duty_grid(d);
info.phi = phi;
info.phi_best = pmin;
info.phi_no_timing = score(rec, Y, nus, zeros(size(stim_idx(:))));
info.duty_grid = duty_grid;
end

function p = score(rec, Y, nus, tau)
p = 0;
for k = 1:numel(nus)
    [~, pk] = kt_solve(rec, Y(:, k), tau, nus(k));
    p = p + pk;
end
end
