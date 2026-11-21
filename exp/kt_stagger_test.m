function out = kt_stagger_test(cfg)
%KT_STAGGER_TEST R2 review: are the voltages of one injection recorded together?
%   out = kt_stagger_test(cfg)
%   The schedule model assumes that tau depends on the injection only
%   ('per_stim'). This test spreads the voltages of an injection over the
%   fraction kappa of its slot ('mixed' timing model; kappa = 0 is per_stim,
%   kappa = 1 is per_meas) and refits duty, shift and direction for every
%   kappa on the whole record and in a two-fold cross-validation (fit on one
%   half, minimum MAP objective evaluated on the other half).
%   Components: the unaliased, isolated lines of the record (harmonics of f0,
%   cardiac line, extra lines), as in kt_exp_decim at the full rate.
%   Logged (rows): phi (whole record), duty, shift, heldout_phi (per fold),
%   all in chi-square units; compare each kappa with kappa = 0.
%
% Author: Hoang Nhut Huynh

run = kt_run_begin(cfg, 'R2stagger');
try
    kt_init_eidors(cfg);
    Mdl = kt_build_models(cfg, false);
    I = Mdl.inv;
    fs = cfg.real.fs;
    [v, pairs] = kt_read_goe_raw(cfg.real.files, cfg.acq.n_elec);
    mrow = kt_match_rows(I.fmdl, pairs);
    I.J = I.J(mrow, :); I.stim_idx = I.stim_idx(mrow); I.k_idx = pairs(:, 3);
    dv = v ./ mean(v, 2) - 1;
    F = size(dv, 2);

    % components (same selection as kt_exp_decim)
    g = sum(dv, 1);
    fgrid = 0:fs/F/4:fs/2;
    G = abs(kt_spectrum(g, fs, fgrid, cfg.rec.window));
    band = fgrid >= cfg.real.f0_search(1) & fgrid <= cfg.real.f0_search(2);
    fb = fgrid(band); [~, ib] = max(G(band)); f0 = fb(ib);
    nu = f0 * (1:cfg.real.n_harm);
    if ~isempty(cfg.real.cardiac_search)
        bandc = fgrid >= cfg.real.cardiac_search(1) & fgrid <= cfg.real.cardiac_search(2);
        fbh = abs(mod(nu + fs/2, fs) - fs/2);
        for k = 1:numel(fbh), bandc = bandc & abs(fgrid - fbh(k)) > 2*fs/F; end
        fc = fgrid(bandc); [~, ic] = max(G(bandc)); nu(end+1) = fc(ic);
    end
    for b = 1:numel(cfg.real.extra_search)
        fb_h = abs(mod(nu + fs/2, fs) - fs/2);
        bandx = fgrid >= cfg.real.extra_search{b}(1) & fgrid <= cfg.real.extra_search{b}(2);
        for k = 1:numel(fb_h), bandx = bandx & abs(fgrid - fb_h(k)) > 2*fs/F; end
        fx = fgrid(bandx); [~, ix] = max(G(bandx)); nu(end+1) = fx(ix); %#ok<AGROW>
    end
    C = kt_alias_map(nu, fs, cfg.rec.tol_factor*fs/F);
    ok = [C.valid] & arrayfun(@(c) sum([C.group] == c.group) == 1, C);
    C = C(ok);
    w = kt_noise_floor(dv, fs, [[C.bin], cfg.real.extra_excl], cfg.rec.window, 4*fs/F);
    prior = cfg.rec.prior;
    if strcmp(prior, 'laplace'), prior = 'external'; end
    rec = kt_recon_matrix(I.J, w, prior, cfg.rec.alpha, Mdl.RtR_laplace);
    fprintf('stagger test: %d components, f0 = %.3f Hz\n', numel(C), f0);

    kappas = kt_get(cfg, 'stagger.kappa', 0:0.25:1);
    halves = {1:floor(F/2), floor(F/2)+1:F};
    Yall = kt_spectrum(dv, fs, [C.bin], cfg.rec.window);
    Yh = cellfun(@(h) kt_spectrum(dv(:, h), fs, [C.bin], cfg.rec.window), halves, 'UniformOutput', false);
    rows = cell(0, 6);
    for kap = kappas
        P0 = struct('n_stim', I.n_stim, 'mode', 'mixed', 'kappa', kap, 'duty', 1, 'shift', 0, 'dir', 1);
        [P, info] = kt_estimate_timing_map(Yall, [C.signed], I.stim_idx, I.k_idx, fs, P0, cfg.rec.duty_grid, rec);
        rows = [rows; {'kappa', kap, 'timing', 'mixed', 'phi', info.phi_best; ...
                       'kappa', kap, 'timing', 'mixed', 'phi_no_timing', info.phi_no_timing; ...
                       'kappa', kap, 'timing', 'mixed', 'duty', P.duty; ...
                       'kappa', kap, 'timing', 'mixed', 'shift', P.shift*P.dir}]; %#ok<AGROW>
        for h = 1:2
            Pa = kt_estimate_timing_map(Yh{h}, [C.signed], I.stim_idx, I.k_idx, fs, P0, cfg.rec.duty_grid, rec);
            tau = kt_timing_model(I.stim_idx, I.k_idx, fs, Pa);
            pb = 0;
            for k = 1:numel(C)
                [~, pk] = kt_solve(rec, Yh{3-h}(:, k), tau, C(k).signed);
                pb = pb + pk;
            end
            rows = [rows; {'kappa', kap, sprintf('fold%d', h), 'mixed', 'heldout_phi', pb}]; %#ok<AGROW>
        end
        fprintf('kappa %.2f: phi %.1f (no timing %.1f), duty %.2f shift %d dir %+d\n', ...
                kap, info.phi_best, info.phi_no_timing, P.duty, P.shift, P.dir);
    end
    out.rows = rows;
    kt_run_end(run, rows, sprintf('stagger test, %d kappa values', numel(kappas)));
catch err
    kt_run_end(run, {}, '', ['FAILED: ' err.message]);
    rethrow(err);
end
end
