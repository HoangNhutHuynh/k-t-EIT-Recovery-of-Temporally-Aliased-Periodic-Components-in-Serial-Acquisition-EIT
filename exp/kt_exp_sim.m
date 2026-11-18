function out = kt_exp_sim(cfg)
%KT_EXP_SIM E1: one simulated sequential acquisition, all methods, logged.
%   out = kt_exp_sim(cfg)
%   Claim tested: with the per-measurement acquisition times in the forward
%   model, component images (fundamental and harmonics above the frame-rate
%   Nyquist) match the ideal simultaneous-acquisition reconstruction, while
%   ignoring timing (naive) or correcting only below Nyquist (yerworth) does not.
%
% Author: Hoang Nhut Huynh

run = kt_run_begin(cfg, 'E1sim');
try
    kt_init_eidors(cfg);
    Mdl = kt_build_models(cfg, true);
    out = kt_sim_core(cfg, Mdl);
    if cfg.save_mat, save(fullfile(run.dir, 'results.mat'), 'out', '-v7'); end

    m = out.methods;
    imgs = [{out.Xt}, cellfun(@(k) out.A.X.(k), m, 'UniformOutput', false)];
    lab = cellfun(@(n, f) sprintf('%s %.1f Hz', n, f), out.names, num2cell(out.nu), 'UniformOutput', false);
    fh = kt_plot_grid(Mdl.inv, imgs, [{'truth'}, m], lab, [true, false(1, numel(m))], ...
        sprintf('E1 |X|  fs=%g Hz, duty=%.2f, timing=%s, SNR=%g dB', cfg.acq.fs, cfg.acq.duty, ...
        cfg.rec.timing_source, cfg.sim.snr_db), 1 + find(strcmp(m, 'simref')));
    kt_save_fig(run, fh, 'fig_components');

    fh = figure('Color', 'w');
    for j = 1:2
        met = {'rho_ref', 'rho_true'};
        subplot(1, 2, j);
        bar(metric_matrix(out.rows, out.names, m, met{j}).');
        set(gca, 'XTickLabel', out.names, 'TickLabelInterpreter', 'none');
        ylim([0 1.05]); ylabel(met{j}, 'Interpreter', 'none'); grid on;
        if j == 2, legend(m, 'Interpreter', 'none', 'Location', 'southoutside', 'Orientation', 'horizontal'); end
    end
    kt_save_fig(run, fh, 'fig_rho');

    S = metric_matrix(out.rows, out.names, m, 'rho_ref');
    summary = sprintf('rho_ref [%s] | kt=%s yerworth=%s naive=%s', strjoin(out.names, ' '), ...
        mat2str(S(strcmp(m, 'kt'), :), 3), mat2str(S(strcmp(m, 'yerworth'), :), 3), ...
        mat2str(S(strcmp(m, 'naive'), :), 3));
    kt_run_end(run, out.rows, summary);
catch err
    kt_run_end(run, {}, '', ['FAILED: ' err.message]);
    rethrow(err);
end
end

function S = metric_matrix(rows, names, methods, metric)
S = nan(numel(methods), numel(names));
for i = 1:size(rows, 1)
    if strcmp(rows{i, 5}, metric)
        a = find(strcmp(methods, rows{i, 4}));
        b = find(strcmp(names, rows{i, 3}));
        if ~isempty(a) && ~isempty(b), S(a, b) = rows{i, 6}; end
    end
end
end
