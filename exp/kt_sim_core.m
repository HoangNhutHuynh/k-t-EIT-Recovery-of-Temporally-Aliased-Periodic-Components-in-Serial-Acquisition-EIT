function out = kt_sim_core(cfg, Mdl)
%KT_SIM_CORE One simulated record: acquire, reconstruct with every method, score.
%   out = kt_sim_core(cfg, Mdl)   (no logging, no figures; used by E1-E3)
%   Scores (out.rows, columns {param, pvalue, component, method, metric, value}):
%     *_ref  vs the reconstruction of ideal simultaneous data with the same
%            solver ('simref'): isolates the effect of the timing model
%     *_true vs the phantom on the inverse mesh (includes regularisation blur;
%            the only valid reference when aliases collide, because
%            simultaneous sampling cannot separate them)
%     quad_tv_ref : total-variation distance of the quadrant fractions to simref
%     phi         : MAP objective (data term in chi-square units)
%     nu_used     : signed frequency the method assumed
%
% Author: Hoang Nhut Huynh

fs = cfg.acq.fs;
I = Mdl.inv;
P_true = struct('n_stim', I.n_stim, 'mode', cfg.acq.timing, 'duty', cfg.acq.duty, ...
                'shift', cfg.acq.shift, 'dir', cfg.acq.dir);
tau_true = kt_timing_model(Mdl.fwd.stim_idx, Mdl.fwd.k_idx, fs, P_true);
[dv_seq, dv_sim] = kt_simulate(cfg, Mdl, tau_true);
F = size(dv_seq, 2);

[nu_all, Xt_all, names_all, analyze] = kt_phantom(cfg, I.cent, I.R);
nu = nu_all(analyze);
Xt = Xt_all(:, analyze);
names = names_all(analyze);

opt = struct('window', cfg.rec.window, 'tol', 2*fs/F, 'l_max', cfg.rec.l_max);
opt.nu_max = inf(1, numel(nu)); opt.nu_max(strcmp(names, 'cardiac')) = kt_get(cfg, 'rec.cardiac_numax', inf);
Call = kt_alias_map(nu_all, fs, opt.tol);
w = kt_noise_floor(dv_seq, fs, [Call.bin], cfg.rec.window, 4*fs/F);
prior = cfg.rec.prior;
if strcmp(prior, 'laplace'), prior = 'external'; end
rec = kt_recon_matrix(I.J, w, prior, cfg.rec.alpha, Mdl.RtR_laplace);

rows = cell(0, 6);
switch cfg.rec.timing_source
    case 'true'
        P = P_true;
    case 'assumed'
        P = P_true;
        f = fieldnames(cfg.rec.timing_assumed);
        for k = 1:numel(f), P.(f{k}) = cfg.rec.timing_assumed.(f{k}); end
    case 'estimate'
        f0 = nu(strcmp(names, cfg.rec.timing_ref));
        y0 = kt_spectrum(dv_seq, fs, f0, cfg.rec.window);
        P0 = P_true; P0.mode = cfg.rec.timing_assumed.mode;
        [Pc, ti] = kt_estimate_timing(y0, f0, I.stim_idx, I.k_idx, fs, P0, cfg.rec.duty_grid, w);
        fprintf('timing (coherence): duty %.2f shift %d dir %+d (true %.2f %d %+d), coherence %.4f (ignored %.4f)\n', ...
            Pc.duty, Pc.shift, Pc.dir, P_true.duty, P_true.shift, P_true.dir, ti.c_best, ti.c_no_timing);
        Ct = kt_alias_map(nu, fs, opt.tol);
        ok = [Ct.valid] & arrayfun(@(c) sum([Ct.group] == c.group) == 1, Ct);   % skip colliding bins
        Yt = kt_spectrum(dv_seq, fs, [Ct(ok).bin], cfg.rec.window);
        [Pm, tm] = kt_estimate_timing_map(Yt, [Ct(ok).signed], I.stim_idx, I.k_idx, fs, P0, cfg.rec.duty_grid, rec);
        fprintf('timing (MAP, %d components): duty %.2f shift %d dir %+d, dPhi vs no timing %.1f\n', ...
            nnz(ok), Pm.duty, Pm.shift, Pm.dir, tm.phi_no_timing - tm.phi_best);
        rows = [rows; {'', NaN, 'timing', 'coherence', 'duty_err', Pc.duty - P_true.duty; ...
                       '', NaN, 'timing', 'coherence', 'shift_ok', double(Pc.shift == P_true.shift && Pc.dir == P_true.dir)}];
        switch cfg.rec.timing_method
            case 'map',       P = Pm;
            case 'coherence', P = Pc;
            otherwise, error('kt_sim_core: unknown timing_method "%s"', cfg.rec.timing_method);
        end
        rows = [rows; {'', NaN, 'timing', 'estimate', 'duty_err', P.duty - P_true.duty; ...
                       '', NaN, 'timing', 'estimate', 'shift_ok', double(P.shift == P_true.shift && P.dir == P_true.dir); ...
                       '', NaN, 'timing', 'estimate', 'coherence', ti.c_best; ...
                       '', NaN, 'timing', 'estimate', 'coherence_no_timing', ti.c_no_timing}];
    otherwise
        error('kt_sim_core: unknown timing_source "%s"', cfg.rec.timing_source);
end
tau = kt_timing_model(I.stim_idx, I.k_idx, fs, P);

A = kt_analyze(dv_seq, fs, nu, tau, rec, cfg.rec.methods, opt);
R0 = kt_analyze(dv_sim, fs, nu, tau, rec, {'naive'}, opt);
A.X.simref = R0.X.naive;
A.Phi.simref = R0.Phi.naive;
A.nu_used.simref = nan(size(nu));
methods = [cfg.rec.methods, {'simref'}];
if cfg.rec.try_eidors_ftc
    Xf = kt_try_ftc(dv_seq, fs, nu, tau, rec, opt);
    if ~isempty(Xf)
        A.X.ftc = Xf; A.Phi.ftc = nan(size(nu)); A.nu_used.ftc = nan(size(nu));
        methods = [methods, {'ftc'}];
    end
end

for c = 1:numel(nu)
    if ~A.C(c).valid
        fprintf('component %s (%.2f Hz) folds onto DC/Nyquist at fs = %g Hz: skipped\n', names{c}, nu(c), fs);
        continue
    end
    qref = kt_quadrants(A.X.simref(:, c), I.cent, I.area);
    for k = 1:numel(methods)
        m = methods{k};
        x = A.X.(m)(:, c);
        sr = kt_compare(x, A.X.simref(:, c), I.area);
        st = kt_compare(x, Xt(:, c), I.area);
        q = kt_quadrants(x, I.cent, I.area);
        rows = [rows; {'', NaN, names{c}, m, 'rho_ref', sr.rho; ...
                       '', NaN, names{c}, m, 'amp_ref', sr.amp; ...
                       '', NaN, names{c}, m, 'relerr_ref', sr.relerr; ...
                       '', NaN, names{c}, m, 'rho_true', st.rho; ...
                       '', NaN, names{c}, m, 'relerr_true', st.relerr; ...
                       '', NaN, names{c}, m, 'quad_tv_ref', 0.5*sum(abs(q - qref)); ...
                       '', NaN, names{c}, m, 'phi', A.Phi.(m)(c); ...
                       '', NaN, names{c}, m, 'nu_used', A.nu_used.(m)(c)}]; %#ok<AGROW>
    end
end

out.rows = rows;
out.A = A;
out.methods = methods;
out.names = names;
out.nu = nu;
out.Xt = Xt;
out.rec = rec;
out.tau = tau;
out.P = P;
end
