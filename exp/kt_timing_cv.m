function rows = kt_timing_cv(dv, fs, C, ok, I, rec, P0, cfg, do_free)
%KT_TIMING_CV Two-fold check that the estimated schedule is not over-fitted.
%   rows = kt_timing_cv(dv, fs, C, ok, I, rec, P0, cfg, do_free)
%   The record is split into halves A and B. The schedule is estimated on one
%   half (equal-slot MAP, and optionally free per-injection timing) and the
%   objective is evaluated on the other half:
%     dphi_heldout_uniform    = Phi_B(no timing) - Phi_B(equal-slot from A)  (> 0: timing generalises)
%     dphi_heldout_free       = Phi_B(equal-slot from A) - Phi_B(free from A) (> 0: free timing generalises)
%     dphi_insample_free      = same on the half used for fitting (optimistic)
%   Phi in chi-square units (weights scaled to the half length).
%
% Author: Hoang Nhut Huynh

F = size(dv, 2);
h = floor(F/2);
part = {1:h, h+1:2*h};
s = h/F;
recs = rec; recs.w = rec.w*s; recs.H = rec.H*s; recs.lam2 = rec.lam2*s;
nus = [C(ok).signed];
bins = [C(ok).bin];
rows = cell(0, 6);
for a = 1:2
    b = 3 - a;
    Ya = kt_spectrum(dv(:, part{a}), fs, bins, cfg.rec.window);
    Yb = kt_spectrum(dv(:, part{b}), fs, bins, cfg.rec.window);
    Pa = kt_estimate_timing_map(Ya, nus, I.stim_idx, I.k_idx, fs, P0, cfg.rec.duty_grid, recs);
    tauA = kt_timing_model(I.stim_idx, I.k_idx, fs, Pa);
    phi0 = score(recs, Yb, nus, zeros(size(tauA)));
    phiU = score(recs, Yb, nus, tauA);
    m = sprintf('fold%d', a);
    rows = [rows; {'', NaN, 'timing_cv', m, 'duty_est', Pa.duty; ...
                   '', NaN, 'timing_cv', m, 'dphi_heldout_uniform', phi0 - phiU}]; %#ok<AGROW>
    fprintf('timing CV fold %d: equal-slot duty %.2f from half %d, held-out dPhi vs no timing %.1f\n', ...
        a, Pa.duty, a, phi0 - phiU);
    if do_free
        tauF = kt_estimate_timing_free(Ya, nus, I.stim_idx, fs, recs, tauA, 4, 60);
        phiF = score(recs, Yb, nus, tauF);
        inF = score(recs, Ya, nus, tauF);
        inU = score(recs, Ya, nus, tauA);
        rows = [rows; {'', NaN, 'timing_cv', m, 'dphi_heldout_free', phiU - phiF; ...
                       '', NaN, 'timing_cv', m, 'dphi_insample_free', inU - inF}]; %#ok<AGROW>
        fprintf('timing CV fold %d: free vs equal-slot, in-sample %.1f, held-out %.1f\n', a, inU - inF, phiU - phiF);
    end
end
end

function p = score(rec, Y, nus, tau)
p = 0;
for k = 1:numel(nus)
    [~, pk] = kt_solve(rec, Y(:, k), tau, nus(k));
    p = p + pk;
end
end
