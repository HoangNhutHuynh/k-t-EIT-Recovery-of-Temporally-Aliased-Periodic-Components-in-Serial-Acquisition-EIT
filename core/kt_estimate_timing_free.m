function [tau, info] = kt_estimate_timing_free(Y, nus, stim_idx, fs, rec, tau0, n_sweep, n_grid)
%KT_ESTIMATE_TIMING_FREE One acquisition time per injection, by minimum MAP objective.
%   [tau, info] = kt_estimate_timing_free(Y, nus, stim_idx, fs, rec, tau0, n_sweep, n_grid)
%   Drops the equal-slot assumption of kt_timing_model: every injection s has
%   its own offset t_s in [0, 1.5/fs). Coordinate descent: for each injection,
%   grid search over n_grid values with the others fixed; n_sweep passes.
%   The first injection present is the gauge (t = 0); a common offset only
%   changes the global phase. Score as kt_estimate_timing_map (sum of Phi_k).
%   tau0 : M x 1 start (e.g. from kt_estimate_timing_map)
%   info.t_s (offset per injection, in frame periods), info.phi per sweep.
%
% Author: Hoang Nhut Huynh

inj = unique(stim_idx(:)).';
t = zeros(1, max(inj));
for s = inj
    t(s) = mean(tau0(stim_idx == s));
end
t = t - t(inj(1));
grid = (0:n_grid-1) / n_grid * 1.5 / fs;
phi_sweep = zeros(1, n_sweep);
for it = 1:n_sweep
    for s = inj(2:end)
        best = inf;
        for g = grid
            t(s) = g;
            p = score(rec, Y, nus, t(stim_idx(:)).');
            if p < best, best = p; gbest = g; end
        end
        t(s) = gbest;
    end
    phi_sweep(it) = best;
    fprintf('  free timing sweep %d: Phi = %.1f\n', it, best);
end
tau = t(stim_idx(:)).';
info.t_s = t * fs;
info.phi = phi_sweep;
end

function p = score(rec, Y, nus, tau)
p = 0;
for k = 1:numel(nus)
    [~, pk] = kt_solve(rec, Y(:, k), tau, nus(k));
    p = p + pk;
end
end
