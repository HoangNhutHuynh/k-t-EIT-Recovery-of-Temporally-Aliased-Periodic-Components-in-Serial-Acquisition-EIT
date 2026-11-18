function out = kt_exp_collision(cfg, duties)
%KT_EXP_COLLISION E3: two aliases in one bin (f0 = fs/4: 3*f0 folds onto -f0).
%   out = kt_exp_collision(cfg, duties)
%   Claims tested:
%   (a) simultaneous sampling cannot separate the aliases (principal angle 0),
%       sequential sampling can, and separability grows with the spread of
%       acquisition times (duty) and shrinks with the spatial rank kept;
%   (b) the joint solve ('kt') recovers the colliding harmonic, the
%       one-at-a-time solve ('kt_single') does not. Scored against the truth.
%   cfg.sim.f0 is set to fs/4 here; components = fundamental (lungs),
%   3rd harmonic (right lung only), cardiac, drift. The 2nd harmonic is
%   exactly at Nyquist and is left out.
%
% Author: Hoang Nhut Huynh

if nargin < 2, duties = [0.25 0.5 0.75 0.9 1.0]; end
fs = cfg.acq.fs;
cfg.sim.f0 = fs/4;
keep = ismember({cfg.sim.components.name}, {'hfov_f0', 'hfov_h3', 'cardiac', 'drift'});
cfg.sim.components = cfg.sim.components(keep);
cfg.sim.components(strcmp({cfg.sim.components.name}, 'hfov_h3')).amp = -0.02;
cfg.rec.methods = {'naive', 'yerworth', 'kt_single', 'kt'};
cfg.rec.timing_source = 'true';
cfg.rec.try_eidors_ftc = false;

run = kt_run_begin(cfg, 'E3collision');
try
    kt_init_eidors(cfg);
    Mdl = kt_build_models(cfg, true);
    rows = cell(0, 6);
    for i = 1:numel(duties)
        c_i = cfg;
        c_i.acq.duty = duties(i);
        o = kt_sim_core(c_i, Mdl);
        r = o.rows; r(:, 1) = {'acq.duty'}; r(:, 2) = {duties(i)};
        rows = [rows; r]; %#ok<AGROW>
        if abs(duties(i) - cfg.acq.duty) < 1e-9, out = o; end
    end
    if ~exist('out', 'var'), out = o; end

    % separability of the pair (+f0, -3 f0): dnu = -fs
    ranks = 5:5:min(200, size(Mdl.inv.J, 1));
    modes = {'per_stim', 'per_meas'};
    fh = figure('Color', 'w', 'Position', [50 50 900 330]);
    for a = 1:2
        subplot(1, 2, a); hold on;
        for i = 1:numel(duties)
            P = struct('n_stim', Mdl.inv.n_stim, 'mode', modes{a}, 'duty', duties(i), ...
                       'shift', cfg.acq.shift, 'dir', cfg.acq.dir);
            tau = kt_timing_model(Mdl.inv.stim_idx, Mdl.inv.k_idx, fs, P);
            s = kt_separability(Mdl.inv.J, out.rec.w, tau, -fs, ranks);
            plot(s.ranks, s.min_angle_deg, '-', 'DisplayName', sprintf('duty %.2f', duties(i)));
            for k = 1:numel(s.ranks)
                rows = [rows; {sprintf('sep_%s_duty', modes{a}), duties(i), 'pair_f0_3f0', ...
                               sprintf('rank%d', s.ranks(k)), 'min_angle_deg', s.min_angle_deg(k)}]; %#ok<AGROW>
            end
        end
        xlabel('rank r of W^{1/2}J kept'); ylabel('smallest principal angle (deg)');
        title(['separability, ' modes{a}], 'Interpreter', 'none'); grid on;
        legend('Location', 'northeast');
    end
    kt_save_fig(run, fh, 'fig_separability');

    m = out.methods;
    imgs = [{out.Xt}, cellfun(@(k) out.A.X.(k), m, 'UniformOutput', false)];
    lab = cellfun(@(n, f) sprintf('%s %.1f Hz', n, f), out.names, num2cell(out.nu), 'UniformOutput', false);
    fh = kt_plot_grid(Mdl.inv, imgs, [{'truth'}, m], lab, [true, false(1, numel(m))], ...
        sprintf('E3 collision f0 = %g Hz, fs = %g Hz, duty = %.2f', cfg.sim.f0, fs, out.P.duty), 1 + find(strcmp(m, 'kt')));
    kt_save_fig(run, fh, 'fig_components');
    if cfg.save_mat, save(fullfile(run.dir, 'results.mat'), 'rows', 'out', '-v7'); end

    pick = @(meth) cellfun(@(x) x, rows(strcmp(rows(:, 3), 'hfov_h3') & strcmp(rows(:, 4), meth) & ...
        strcmp(rows(:, 5), 'rho_true'), 6)).';
    summary = sprintf('h3 rho_true vs duty %s: kt=%s kt_single=%s', mat2str(duties, 3), ...
        mat2str(pick('kt'), 3), mat2str(pick('kt_single'), 3));
    kt_run_end(run, rows, summary);
catch err
    kt_run_end(run, {}, '', ['FAILED: ' err.message]);
    rethrow(err);
end
end
