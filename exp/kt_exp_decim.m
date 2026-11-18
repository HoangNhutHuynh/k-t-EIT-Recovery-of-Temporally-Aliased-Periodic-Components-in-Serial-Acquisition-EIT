function out = kt_exp_decim(cfg)
%KT_EXP_DECIM E7: frame-rate reduction on real data, scored against a measured reference.
%   out = kt_exp_decim(cfg)
%   Keeping every D-th frame (D = cfg.decim.factor) lowers the frame rate to
%   fs/D while every measurement keeps its acquisition offset tau inside the
%   frame. Harmonics that are below Nyquist at fs but above it at fs/D become
%   aliased. Their images reconstructed at the full rate (no aliasing; k-t and
%   yerworth coincide there) are the reference, so every method can be scored
%   on real data with rho_ref / amp_ref / relerr_ref, as in the simulations.
%   All D frame offsets are used (D disjoint sub-records); metrics are logged
%   per offset (param 'offset').
%   cfg.decim.timing: 'decimated' re-estimates the schedule (equal slots, MAP)
%   from the low-rate data only; 'full' uses the schedule from the full rate.
%   Components: harmonics of f0 that are unaliased at fs (plus the cardiac line
%   when cfg.real.cardiac_search is set).
%
% Author: Hoang Nhut Huynh

run = kt_run_begin(cfg, 'E7decim');
try
    kt_init_eidors(cfg);
    Mdl = kt_build_models(cfg, false);
    I = Mdl.inv;
    fs = cfg.real.fs;
    [v, pairs, info] = kt_read_goe_raw(cfg.real.files, cfg.acq.n_elec);
    mrow = kt_match_rows(I.fmdl, pairs);
    I.J = I.J(mrow, :); I.vh = I.vh(mrow); I.stim_idx = I.stim_idx(mrow); I.k_idx = pairs(:, 3);
    dv = v ./ mean(v, 2) - 1;
    F = size(dv, 2);
    fprintf('%s: %d frames at %g Hz, %d measurements\n', info.format, F, fs, size(dv, 1));

    % --- components unaliased at the full rate
    g = sum(dv, 1);
    fgrid = 0:fs/F/4:fs/2;
    G = abs(kt_spectrum(g, fs, fgrid, cfg.rec.window));
    band = fgrid >= cfg.real.f0_search(1) & fgrid <= cfg.real.f0_search(2);
    fb = fgrid(band); [~, ib] = max(G(band)); f0 = fb(ib);
    nu = f0 * (1:cfg.real.n_harm);
    names = arrayfun(@(h) sprintf('h%d', h), 1:numel(nu), 'UniformOutput', false);
    if ~isempty(cfg.real.cardiac_search)
        bandc = fgrid >= cfg.real.cardiac_search(1) & fgrid <= cfg.real.cardiac_search(2);
        fbh = abs(mod(nu + fs/2, fs) - fs/2);
        for k = 1:numel(fbh), bandc = bandc & abs(fgrid - fbh(k)) > 2*fs/F; end
        fc = fgrid(bandc); [~, ic] = max(G(bandc));
        nu(end+1) = fc(ic); names{end+1} = 'cardiac';
    end
    for b = 1:numel(cfg.real.extra_search)    % extra unaliased lines (e.g. slow respiration)
        fb_h = abs(mod(nu + fs/2, fs) - fs/2);
        bandx = fgrid >= cfg.real.extra_search{b}(1) & fgrid <= cfg.real.extra_search{b}(2);
        for k = 1:numel(fb_h), bandx = bandx & abs(fgrid - fb_h(k)) > 2*fs/F; end
        fx = fgrid(bandx); [~, ix] = max(G(bandx));
        nu(end+1) = fx(ix); names{end+1} = sprintf('line%d', b); %#ok<AGROW>
        fprintf('extra line %d at %.3f Hz\n', b, fx(ix));
    end
    opt = struct('window', cfg.rec.window, 'tol', cfg.rec.tol_factor*fs/F, 'l_max', cfg.rec.l_max);
    C = kt_alias_map(nu, fs, opt.tol);
    keep = [C.valid] & [C.ell] == 0 & arrayfun(@(c) sum([C.group] == c.group) == 1, C);
    nu = nu(keep); names = names(keep); C = C(keep);
    fprintf('f0 = %.3f Hz; reference components (unaliased at %g Hz): %s\n', f0, fs, strjoin(names, ' '));

    % --- reference at the full rate
    w = kt_noise_floor(dv, fs, [[C.bin], cfg.real.extra_excl], cfg.rec.window, 4*fs/F);
    prior = cfg.rec.prior;
    if strcmp(prior, 'laplace'), prior = 'external'; end
    rec = kt_recon_matrix(I.J, w, prior, cfg.rec.alpha, Mdl.RtR_laplace);
    P0 = struct('n_stim', I.n_stim, 'mode', cfg.rec.timing_assumed.mode, 'duty', 1, 'shift', 0, 'dir', 1);
    Yf = kt_spectrum(dv, fs, [C.bin], cfg.rec.window);
    Pf = kt_estimate_timing_map(Yf, [C.signed], I.stim_idx, I.k_idx, fs, P0, cfg.rec.duty_grid, rec);
    tau_f = kt_timing_model(I.stim_idx, I.k_idx, fs, Pf);
    Aref = kt_analyze(dv, fs, nu, tau_f, rec, {'kt'}, opt);
    Xref = Aref.X.kt;
    fprintf('full-rate schedule: duty %.2f shift %d dir %+d\n', Pf.duty, Pf.shift, Pf.dir);

    % --- reduced frame rate
    D = cfg.decim.factor;
    fsd = fs / D;
    nud = nu * (1 + kt_get(cfg, 'decim.f0_err', 0));   % frequencies used at the reduced rate (sensitivity to f0 error)
    if nud(1) ~= nu(1), fprintf('reduced-rate reconstruction with f0 error %+.2f %%\n', 100*kt_get(cfg, 'decim.f0_err', 0)); end
    rows = {'', NaN, 'timing', 'full_rate', 'duty', Pf.duty; '', NaN, 'signal', 'global', 'f0', f0};
    methods = cfg.rec.methods;
    Cd = kt_alias_map(nud, fsd, cfg.rec.tol_factor*fsd/floor(F/D));
    for c = 1:numel(C)
        fprintf('  %-8s %6.2f Hz -> bin %5.2f Hz at %g Hz (ell %d)%s\n', names{c}, nu(c), Cd(c).bin, fsd, Cd(c).ell, ...
            ifelse(Cd(c).valid, '', '  [DC/Nyquist: skipped]'));
        rows = [rows; {'', NaN, names{c}, 'map', 'ell_decimated', Cd(c).ell}]; %#ok<AGROW>
    end
    for o = 0:D-1
        dvd = dv(:, 1+o:D:end);
        Fd = size(dvd, 2);
        optd = opt; optd.tol = cfg.rec.tol_factor*fsd/Fd;
        Cdo = kt_alias_map(nud, fsd, optd.tol);
        okd = [Cdo.valid] & arrayfun(@(c) sum([Cdo.group] == c.group) == 1, Cdo);
        switch cfg.decim.timing
            case 'full'
                tau_d = tau_f;
            case 'decimated'
                Yd = kt_spectrum(dvd, fsd, [Cdo(okd).bin], cfg.rec.window);
                Pd = kt_estimate_timing_map(Yd, [Cdo(okd).signed], I.stim_idx, I.k_idx, fsd, P0, cfg.rec.duty_grid, rec);
                tau_d = kt_timing_model(I.stim_idx, I.k_idx, fsd, Pd);
                fprintf('offset %d: schedule at %g Hz: duty %.2f (expected %.2f) shift %d dir %+d\n', ...
                    o, fsd, Pd.duty, Pf.duty/D, Pd.shift, Pd.dir);
                rows = [rows; {'offset', o, 'timing', 'decimated', 'duty_x_factor', Pd.duty*D}]; %#ok<AGROW>
            otherwise
                error('kt_exp_decim: unknown decim.timing "%s"', cfg.decim.timing);
        end
        Ad = kt_analyze(dvd, fsd, nud, tau_d, rec, methods, optd);
        mlist = methods;
        if cfg.rec.try_eidors_ftc
            Xf = kt_try_ftc(dvd, fsd, nud, tau_d, rec, optd, kt_fill_reciprocal(dvd, pairs, Mdl.inv.fmdl), mrow, kt_get(cfg, 'rec.ftc_method', 'FFT'));
            if ~isempty(Xf), Ad.X.ftc = Xf; mlist = [mlist, {'ftc'}]; end %#ok<AGROW>
        end
        for k = 1:numel(cfg.rec.eidors_temporal)
            mt = cfg.rec.eidors_temporal{k};
            if strcmp(mt, 'time_prior')
                Xt = kt_time_prior(rec, dvd, fsd, [Cdo.bin], cfg.rec.window, 2, 0.8);
            else
                Xt = kt_eidors_temporal(Mdl.imdl, Mdl.inv, mrow, dvd, fsd, [Cdo.bin], cfg.rec.window, mt, cfg.rec.alpha);
            end
            if ~isempty(Xt), Ad.X.(mt) = Xt; mlist = [mlist, {mt}]; end %#ok<AGROW>
        end
        if o == 0, A0 = Ad; m0 = mlist; end
        for c = 1:numel(C)
            if ~Cdo(c).valid, continue; end
            for k = 1:numel(mlist)
                m = mlist{k};
                s = kt_compare(Ad.X.(m)(:, c), Xref(:, c), I.area);
                rows = [rows; {'offset', o, names{c}, m, 'rho_ref', s.rho; 'offset', o, names{c}, m, 'amp_ref', s.amp; ...
                               'offset', o, names{c}, m, 'relerr_ref', s.relerr}]; %#ok<AGROW>
            end
            if isfield(Ad.nu_used, 'kt_auto')
                rows = [rows; {'offset', o, names{c}, 'kt_auto', 'kt_auto_ok', ...
                    double(abs(Ad.nu_used.kt_auto(c) - Cdo(c).signed) < 1e-6)}]; %#ok<AGROW>
            end
        end
    end

    % --- summary: aliased components, mean over offsets
    al = find([Cd.ell] ~= 0 & [Cd.valid]);

    % --- figures: reference (full rate) vs every method at the reduced rate (offset 0)
    order = {'naive', 'yerworth', 'ftc', 'time_prior', 'kt'};
    label = {'Naive', 'Phase corr.', 'FTC', 'Temporal prior', 'k-t (proposed)'};
    keepm = ismember(order, m0);
    show = order(keepm); lab_m = label(keepm);
    if ~isempty(al)
        imgs = [{Xref(:, al)}, cellfun(@(k) A0.X.(k)(:, al), show, 'UniformOutput', false)];
        lab = arrayfun(@(c) sprintf('%s %.1f Hz', names{c}, nu(c)), al, 'UniformOutput', false);
        cols = [{sprintf('Reference (%g fps)', fs)}, lab_m];
        fh = kt_plot_grid(I, imgs, cols, lab, false(1, numel(show) + 1), sprintf('E7 %g -> %g fps', fs, fsd), 1, ...
                          struct('mode', 'aligned'));
        kt_save_fig(run, fh, 'fig_decim_aligned');
        fh = kt_plot_grid(I, imgs, cols, lab, false(1, numel(show) + 1), sprintf('E7 %g -> %g fps', fs, fsd), 1, ...
                          struct('mode', 'abs'));
        kt_save_fig(run, fh, 'fig_decim_abs');
        fh = kt_plot_grid(I, imgs, cols, lab, false(1, numel(show) + 1), sprintf('E7 %g -> %g fps', fs, fsd), 1, ...
                          struct('mode', 'error'));
        kt_save_fig(run, fh, 'fig_decim_error');
        sel = find(ismember(names(al), cfg.decim.fig_rows));
        if ~isempty(sel)                 % main-text figure: selected rows only
            imgs_s = cellfun(@(x) x(:, sel), imgs, 'UniformOutput', false);
            fh = kt_plot_grid(I, imgs_s, cols, lab(sel), false(1, numel(show) + 1), sprintf('E7 %g -> %g fps', fs, fsd), 1, ...
                              struct('mode', 'aligned'));
            kt_save_fig(run, fh, 'fig_decim_main');
            fh = kt_plot_grid(I, imgs_s, cols, lab(sel), false(1, numel(show) + 1), sprintf('E7 %g -> %g fps', fs, fsd), 1, ...
                              struct('mode', 'error'));
            kt_save_fig(run, fh, 'fig_decim_error_main');
        end
    end
    parts = {};
    for c = al
        for m = {'yerworth', 'kt'}
            sel = strcmp(rows(:, 3), names{c}) & strcmp(rows(:, 4), m{1}) & strcmp(rows(:, 5), 'rho_ref');
            parts{end+1} = sprintf('%s %s %.3f', names{c}, m{1}, mean(cell2mat(rows(sel, 6)))); %#ok<AGROW>
        end
    end
    summary = sprintf('fs %g -> %g Hz, timing %s | rho_ref (aliased): %s', fs, fsd, cfg.decim.timing, strjoin(parts, '; '));
    out.rows = rows;
    if cfg.save_mat, save(fullfile(run.dir, 'results.mat'), 'out', 'Xref', '-v7'); end
    kt_run_end(run, rows, summary);
catch err
    kt_run_end(run, {}, '', ['FAILED: ' err.message]);
    rethrow(err);
end
end

function s = ifelse(c, a, b)
if c, s = a; else, s = b; end
end
