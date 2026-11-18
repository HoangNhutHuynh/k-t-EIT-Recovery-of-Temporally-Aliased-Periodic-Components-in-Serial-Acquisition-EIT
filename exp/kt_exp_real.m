function out = kt_exp_real(cfg)
%KT_EXP_REAL E4: HFOV recording (chb_hfov_porcine), no ground truth, logged.
%   out = kt_exp_real(cfg)
%   Steps: read voltages -> normalised differences -> HFOV fundamental f0 from
%   the global signal -> acquisition schedule from the coherence at f0 ->
%   component images of f0..n_harm*f0 for every method -> internal checks:
%     dPhi_vs_yerworth : Phi(yerworth) - Phi(method), chi-square units. For
%                        one-frequency models the normal matrix is the same for
%                        every assumed frequency, so this is twice the log
%                        evidence ratio (> 0: the method explains the data better)
%     corr_f0, roi_frac: harmonic image vs the method's own fundamental image
%                        (harmonics of lung oscillation should sit in the lungs)
%     kt_auto_ok       : the alias order chosen by evidence equals the order
%                        expected from h*f0
%   All block metrics: mean and percentile bootstrap CI over non-overlapping
%   blocks of cfg.real.block_sec.
%
% Author: Hoang Nhut Huynh

run = kt_run_begin(cfg, 'E4real');
try
    kt_init_eidors(cfg);
    Mdl = kt_build_models(cfg, false);
    I = Mdl.inv;
    fs = cfg.real.fs;
    if strcmp(cfg.real.reader, 'goe_raw')
        [v, pairs, info] = kt_read_goe_raw(cfg.real.files, cfg.acq.n_elec);
        mrow = kt_match_rows(I.fmdl, pairs);        % model row of every stored measurement
        I.J = I.J(mrow, :); I.vh = I.vh(mrow); I.stim_idx = I.stim_idx(mrow);
        I.k_idx = pairs(:, 3);                      % position inside the injection (device order)
        fprintf('model restricted to the %d stored measurements\n', numel(mrow));
    else
        [v, info] = kt_load_real(cfg.real, cfg.acq.n_elec);
        mrow = (1:size(v, 1)).'; pairs = [];
    end
    if isfield(I, 'lung_poly')                      % model check figure (always saved, small)
        fh = figure('Color', 'w');
        patch('Faces', I.elems, 'Vertices', I.nodes, 'FaceColor', 'none', 'EdgeColor', [0.8 0.8 0.8]); hold on;
        for k = 1:numel(I.lung_poly), plot(I.lung_poly{k}([1:end 1], 1), I.lung_poly{k}([1:end 1], 2), 'b-', 'LineWidth', 1.5); end
        for e = 1:numel(I.fmdl.electrode)
            p = mean(I.fmdl.nodes(I.fmdl.electrode(e).nodes, 1:2), 1);
            text(p(1)*1.08, p(2)*1.08, sprintf('%d', e), 'HorizontalAlignment', 'center');
        end
        axis equal off; title(sprintf('%s (electrode 1 = sternum)', cfg.model.shape), 'Interpreter', 'none');
        run_fig = run; run_fig.cfg.save_figs = true;
        kt_save_fig(run_fig, fh, 'fig_model');
    end
    i0 = max(1, floor(cfg.real.t_range(1)*fs) + 1);
    i1 = min(size(v, 2), floor(cfg.real.t_range(2)*fs));
    v = v(:, i0:i1);
    mv = mean(v, 2);
    if any(abs(mv) < 1e-9)
        error('kt_exp_real: %d channels have zero mean (unstored measurements?); use cfg.real.reader = ''goe_raw''', nnz(abs(mv) < 1e-9));
    end
    dv = v ./ mv - 1;
    F = size(dv, 2);
    fprintf('segment %.1f-%.1f s: %d frames, format %s\n', (i0-1)/fs, i1/fs, F, info.format);

    % --- HFOV fundamental from the global signal
    g = sum(dv, 1);
    fgrid = 0:fs/F/4:fs/2;
    G = abs(kt_spectrum(g, fs, fgrid, cfg.rec.window));
    band = fgrid >= cfg.real.f0_search(1) & fgrid <= cfg.real.f0_search(2);
    fb = fgrid(band); [~, ib] = max(G(band)); f0 = fb(ib);
    nu = f0 * (1:cfg.real.n_harm);
    names = arrayfun(@(h) sprintf('h%d', h), 1:numel(nu), 'UniformOutput', false);
    if ~isempty(cfg.real.cardiac_search)     % unaliased cardiac line: global-spectrum peak away from the harmonic bins
        fb_h = abs(mod(nu + fs/2, fs) - fs/2);
        bandc = fgrid >= cfg.real.cardiac_search(1) & fgrid <= cfg.real.cardiac_search(2);
        for k = 1:numel(fb_h), bandc = bandc & abs(fgrid - fb_h(k)) > 2*fs/F; end
        fc = fgrid(bandc); [~, ic] = max(G(bandc)); fcard = fc(ic);
        nu(end+1) = fcard; names{end+1} = 'cardiac';
        fprintf('cardiac line at %.3f Hz (nearest harmonic bin %.3f Hz)\n', fcard, fb_h(find(abs(fb_h - fcard) == min(abs(fb_h - fcard)), 1)));
    end
    for h = cfg.real.extra_fold_of(:).'      % unaliased line at the folded bin of harmonic h (e.g. cardiac)
        nu(end+1) = abs(mod(h*f0 + fs/2, fs) - fs/2); %#ok<AGROW>
        names{end+1} = sprintf('x_at_h%d', h); %#ok<AGROW>
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
    opt.nu_max = inf(1, numel(nu)); opt.nu_max(strcmp(names, 'cardiac')) = kt_get(cfg, 'rec.cardiac_numax', inf);
    C = kt_alias_map(nu, fs, opt.tol);
    fprintf('f0 = %.3f Hz; harmonics fold to %s Hz (alias orders %s)\n', f0, mat2str([C.bin], 4), mat2str([C.ell]));

    fh = figure('Color', 'w', 'Position', [50 50 900 520]);
    subplot(2, 1, 1); plot((0:F-1)/fs + (i0-1)/fs, g); xlabel('time (s)'); ylabel('sum of dv'); grid on;
    title(sprintf('global signal (%s)', info.format), 'Interpreter', 'none');
    subplot(2, 1, 2); semilogy(fgrid, G); hold on;
    for c = 1:numel(C), xline_(C(c).bin, sprintf('%gf0', c)); end
    xlabel('frequency (Hz, folded)'); ylabel('|DTFT|'); grid on;
    kt_save_fig(run, fh, 'fig_global');

    % --- noise weights, reconstruction matrix
    w = kt_noise_floor(dv, fs, [[C.bin], cfg.real.extra_excl], cfg.rec.window, 4*fs/F);
    prior = cfg.rec.prior;
    if strcmp(prior, 'laplace'), prior = 'external'; end
    rec = kt_recon_matrix(I.J, w, prior, cfg.rec.alpha, Mdl.RtR_laplace);

    % --- acquisition schedule
    rows = cell(0, 6);
    P0 = struct('n_stim', I.n_stim, 'mode', cfg.rec.timing_assumed.mode, 'duty', cfg.rec.timing_assumed.duty, ...
                'shift', cfg.rec.timing_assumed.shift, 'dir', cfg.rec.timing_assumed.dir);
    y0 = kt_spectrum(dv, fs, f0, cfg.rec.window);
    [Pe, ti] = kt_estimate_timing(y0, f0, I.stim_idx, I.k_idx, fs, P0, cfg.rec.duty_grid, w);
    fprintf('timing estimate: duty %.2f shift %d dir %+d, coherence %.4f (timing ignored %.4f)\n', ...
        Pe.duty, Pe.shift, Pe.dir, ti.c_best, ti.c_no_timing);
    rows = [rows; {'', NaN, 'timing', 'estimate', 'duty', Pe.duty; '', NaN, 'timing', 'estimate', 'shift', Pe.shift; ...
                   '', NaN, 'timing', 'estimate', 'dir', Pe.dir; '', NaN, 'timing', 'estimate', 'coherence', ti.c_best; ...
                   '', NaN, 'timing', 'estimate', 'coherence_no_timing', ti.c_no_timing; ...
                   '', NaN, 'signal', 'global', 'f0', f0}];
    ok = [C.valid] & arrayfun(@(c) sum([C.group] == c.group) == 1, C);
    Yt = kt_spectrum(dv, fs, [C(ok).bin], cfg.rec.window);
    [Pm, tm] = kt_estimate_timing_map(Yt, [C(ok).signed], I.stim_idx, I.k_idx, fs, P0, cfg.rec.duty_grid, rec);
    fprintf('timing (MAP, %d harmonics): duty %.2f shift %d dir %+d, dPhi vs no timing %.1f\n', ...
        nnz(ok), Pm.duty, Pm.shift, Pm.dir, tm.phi_no_timing - tm.phi_best);
    rows = [rows; {'', NaN, 'timing', 'map', 'duty', Pm.duty; '', NaN, 'timing', 'map', 'shift', Pm.shift; ...
                   '', NaN, 'timing', 'map', 'dir', Pm.dir; '', NaN, 'timing', 'map', 'dphi_vs_no_timing', tm.phi_no_timing - tm.phi_best}];
    switch cfg.rec.timing_source
        case 'assumed', P = P0;
        otherwise   % 'estimate' (and 'true', which has no meaning here)
            if strcmp(cfg.rec.timing_method, 'coherence'), P = Pe; else, P = Pm; end
    end
    fprintf('reconstruction uses timing_source = %s: duty %.2f shift %d dir %+d\n', ...
        cfg.rec.timing_source, P.duty, P.shift, P.dir);
    tau = kt_timing_model(I.stim_idx, I.k_idx, fs, P);
    if strcmp(cfg.rec.timing_method, 'free') && ~strcmp(cfg.rec.timing_source, 'assumed')
        [tau, tf] = kt_estimate_timing_free(Yt, [C(ok).signed], I.stim_idx, fs, rec, tau, 4, 60);
        fprintf('free timing (frame periods): %s; Phi %.1f (uniform-slot MAP %.1f)\n', ...
            mat2str(round(tf.t_s(unique(I.stim_idx)) * 100) / 100), tf.phi(end), tm.phi_best);
        for s = unique(I.stim_idx(:)).'
            rows = [rows; {'', NaN, 'timing', 'free', sprintf('t_inj%02d', s), tf.t_s(s)}]; %#ok<AGROW>
        end
        rows = [rows; {'', NaN, 'timing', 'free', 'dphi_vs_uniform', tm.phi_best - tf.phi(end)}];
    end
    if cfg.real.cv_timing
        rows = [rows; kt_timing_cv(dv, fs, C, ok, I, rec, P0, cfg, strcmp(cfg.rec.timing_method, 'free'))];
    end
    fh = figure('Color', 'w', 'Position', [50 50 900 300]);
    a = find([1 -1] == Pe.dir);
    subplot(1, 2, 1); plot(ti.duty_grid, squeeze(ti.c(a, Pe.shift + 1, :))); grid on;
    xlabel('duty'); ylabel('coherence at f0'); title(sprintf('dir %+d, shift %d', Pe.dir, Pe.shift));
    subplot(1, 2, 2); plot(0:I.n_stim-1, max(squeeze(ti.c(a, :, :)), [], 2), '-o'); grid on;
    xlabel('shift'); ylabel('max coherence over duty');
    kt_save_fig(run, fh, 'fig_timing');

    % --- full-segment images
    methods = cfg.rec.methods;
    A = kt_analyze(dv, fs, nu, tau, rec, methods, opt);
    if cfg.rec.try_eidors_ftc
        if strcmp(cfg.real.reader, 'goe_raw')      % 208-row model-order frame (reciprocity fills unstored rows)
            Xf = kt_try_ftc(dv, fs, nu, tau, rec, opt, kt_fill_reciprocal(dv, pairs, Mdl.inv.fmdl), mrow);
        else
            Xf = kt_try_ftc(dv, fs, nu, tau, rec, opt);
        end
        if ~isempty(Xf), A.X.ftc = Xf; methods = [methods, {'ftc'}]; end
    end
    for k = 1:numel(cfg.rec.eidors_temporal)
        mt = cfg.rec.eidors_temporal{k};
        if strcmp(mt, 'time_prior')
            Xt = kt_time_prior(rec, dv, fs, [C.bin], cfg.rec.window, 2, 0.8);
        else
            Xt = kt_eidors_temporal(Mdl.imdl, Mdl.inv, mrow, dv, fs, [C.bin], cfg.rec.window, mt, cfg.rec.alpha);
        end
        if ~isempty(Xt), A.X.(mt) = Xt; methods = [methods, {mt}]; end
    end
    if isfield(I, 'lung')
        roi = I.lung;
    else                                           % functional ROI from the k-t fundamental
        x1 = abs(A.X.kt(:, strcmp(names, 'h1')));
        roi = x1 >= 0.3 * max(x1);
    end
    lab = cellfun(@(n, c) sprintf('%s %.1f Hz (bin %.1f)', n, c.nu, c.bin), names, num2cell(C), 'UniformOutput', false);
    imgs = cellfun(@(k) A.X.(k), methods, 'UniformOutput', false);
    fh = kt_plot_grid(I, imgs, methods, lab, false(1, numel(methods)), ...
        sprintf('E4 |X|  f0 = %.2f Hz, duty %.2f shift %d dir %+d', f0, P.duty, P.shift, P.dir), find(strcmp(methods, 'kt')));
    kt_save_fig(run, fh, 'fig_components');
    rows = [rows; score(A, methods, C, I, names, '', fs); kt_nl_index(A, methods, names, C, I, roi, '')];

    % --- blocks
    Fb = round(cfg.real.block_sec * fs);
    nb = floor(F / Fb);
    recb = rec; s = Fb/F;
    recb.w = rec.w*s; recb.H = rec.H*s; recb.lam2 = rec.lam2*s;   % same R, Phi in chi-square units
    brow = cell(0, 6);
    optb = opt; optb.tol = cfg.rec.tol_factor*fs/Fb;   % block resolution decides which bins collide
    for b = 1:nb
        Ab = kt_analyze(dv(:, (b-1)*Fb + (1:Fb)), fs, nu, tau, recb, cfg.rec.methods, optb);
        brow = [brow; score(Ab, cfg.rec.methods, C, I, names, 'blk', fs); kt_nl_index(Ab, cfg.rec.methods, names, C, I, roi, 'blk')]; %#ok<AGROW>
    end
    fprintf('%d blocks of %.1f s\n', nb, cfg.real.block_sec);
    rng(cfg.seed);
    keys = unique(strcat(brow(:, 3), '|', brow(:, 4), '|', brow(:, 5)));
    for k = 1:numel(keys)
        p = strsplit(keys{k}, '|');
        sel = strcmp(brow(:, 3), p{1}) & strcmp(brow(:, 4), p{2}) & strcmp(brow(:, 5), p{3});
        [mu, lo, hi] = kt_boot_ci(cell2mat(brow(sel, 6)), cfg.real.n_boot, 0.95);
        rows = [rows; {'blocks', nb, p{1}, p{2}, [p{3} '_mean'], mu; 'blocks', nb, p{1}, p{2}, [p{3} '_lo'], lo; ...
                       'blocks', nb, p{1}, p{2}, [p{3} '_hi'], hi}]; %#ok<AGROW>
    end

    fh = figure('Color', 'w', 'Position', [50 50 300*numel(nu) 300]);
    mm = cfg.rec.methods(~strcmp(cfg.rec.methods, 'yerworth'));
    for c = 1:numel(nu)
        subplot(1, numel(nu), c); hold on;
        for k = 1:numel(mm)
            g3 = @(suf) get_row(rows, names{c}, mm{k}, ['blk_dPhi_vs_yerworth' suf]);
            mu = g3('_mean'); lo = g3('_lo'); hi = g3('_hi');
            bar(k, mu); plot([k k], [lo hi], 'k-', 'LineWidth', 1.5);
        end
        set(gca, 'XTick', 1:numel(mm), 'XTickLabel', mm, 'TickLabelInterpreter', 'none');
        title(lab{c}, 'Interpreter', 'none'); ylabel('\Delta\Phi vs yerworth (block mean, 95% CI)'); grid on;
    end
    kt_save_fig(run, fh, 'fig_dphi');

    out.A = A; out.C = C; out.f0 = f0; out.P = P; out.timing = ti; out.info = info; out.rows = rows;
    if cfg.save_mat, save(fullfile(run.dir, 'results.mat'), 'out', '-v7'); end
    parts = cell(1, numel(nu));
    for c = 1:numel(nu)
        parts{c} = sprintf('%s:%.1f[%.1f,%.1f]', names{c}, get_row(rows, names{c}, 'kt', 'blk_dPhi_vs_yerworth_mean'), ...
            get_row(rows, names{c}, 'kt', 'blk_dPhi_vs_yerworth_lo'), get_row(rows, names{c}, 'kt', 'blk_dPhi_vs_yerworth_hi'));
    end
    summary = sprintf('f0=%.2f duty=%.2f shift=%d dir=%+d coh=%.3f(no timing %.3f) | dPhi kt-vs-yerworth %s', ...
        f0, P.duty, P.shift, P.dir, ti.c_best, ti.c_no_timing, strjoin(parts, ' '));
    kt_run_end(run, rows, summary);
catch err
    kt_run_end(run, {}, '', ['FAILED: ' err.message]);
    rethrow(err);
end
end

function rows = score(A, methods, C, I, names, pre, fs)
if ~isempty(pre), pre = [pre '_']; end
rows = cell(0, 6);
for c = 1:numel(C)
    if ~C(c).valid, continue; end
    for k = 1:numel(methods)
        m = methods{k};
        x = A.X.(m)(:, c);
        x1 = A.X.(m)(:, 1);
        r = {};
        if isfield(A, 'Phi') && isfield(A.Phi, m) && isfield(A.Phi, 'yerworth')
            r = [r; {'dPhi_vs_yerworth', A.Phi.yerworth(c) - A.Phi.(m)(c)}];
        end
        if c > 1
            r = [r; {'corr_f0', wcorr(abs(x), abs(x1), I.area); 'roi_frac', roi_frac(x, x1, I.area)}];
        end
        if isfield(I, 'lung')
            e = I.area(:) .* abs(x(:)).^2;
            r = [r; {'lung_frac', sum(e(I.lung)) / sum(e)}];
        end
        if strcmp(m, 'kt_auto')
            r = [r; {'kt_auto_ok', double(abs(A.nu_used.kt_auto(c) - C(c).signed) < 1e-6)}];
            r = [r; {'ell_auto', round((A.nu_used.kt_auto(c) - C(c).bin) / fs)}];
        end
        if isempty(pre)
            q = kt_quadrants(x, I.cent, I.area);
            r = [r; {'q1', q(1); 'q2', q(2); 'q3', q(3); 'q4', q(4)}];
        end
        for i = 1:size(r, 1)
            rows = [rows; {'', NaN, names{c}, m, [pre r{i, 1}], r{i, 2}}]; %#ok<AGROW>
        end
    end
end
end

function r = wcorr(a, b, w)
w = w(:) / sum(w);
a = a(:) - sum(w .* a(:)); b = b(:) - sum(w .* b(:));
r = sum(w .* a .* b) / sqrt(sum(w .* a.^2) * sum(w .* b.^2));
end

function f = roi_frac(x, x1, area)
roi = abs(x1) >= 0.3 * max(abs(x1));
e = area(:) .* abs(x(:)).^2;
f = sum(e(roi)) / sum(e);
end

function v = get_row(rows, comp, meth, metric)
i = find(strcmp(rows(:, 3), comp) & strcmp(rows(:, 4), meth) & strcmp(rows(:, 5), metric), 1);
v = NaN;
if ~isempty(i), v = rows{i, 6}; end
end

function xline_(x, s)
yl = ylim;
plot([x x], yl, 'r:');
text(x, yl(2), s, 'Color', 'r', 'VerticalAlignment', 'top');
end
