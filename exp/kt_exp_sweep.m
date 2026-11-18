function rows = kt_exp_sweep(cfg, param, values, metric)
%KT_EXP_SWEEP E2: repeat the simulation over one config parameter, logged as one run.
%   rows = kt_exp_sweep(cfg, param, values, metric)
%   param  : dotted cfg path, e.g. 'sim.snr_db', 'rec.alpha', 'sim.f0',
%            'rec.timing_assumed.duty' (with cfg.rec.timing_source = 'assumed')
%   values : numeric vector
%   metric : plotted metric (default 'rho_ref'); all metrics are logged
%   The noise realisation is the same at every point (paired comparison).
%
% Author: Hoang Nhut Huynh

if nargin < 4, metric = 'rho_ref'; end
if strncmp(param, 'model.', 6) || strcmp(param, 'acq.n_elec')
    error('kt_exp_sweep: parameters that change the mesh are not supported');
end
cfg.tag = sprintf('%s_%s', cfg.tag, strrep(param, '.', '-'));
run = kt_run_begin(cfg, 'E2sweep');
try
    kt_init_eidors(cfg);
    Mdl = kt_build_models(cfg, true);
    rows = cell(0, 6);
    for i = 1:numel(values)
        fprintf('--- %s = %g\n', param, values(i));
        out = kt_sim_core(kt_set(cfg, param, values(i)), Mdl);
        r = out.rows;
        r(:, 1) = {param};
        r(:, 2) = {values(i)};
        rows = [rows; r]; %#ok<AGROW>
    end
    if cfg.save_mat, save(fullfile(run.dir, 'results.mat'), 'rows', '-v7'); end

    names = out.names;
    methods = out.methods;
    fh = figure('Color', 'w', 'Position', [50 50 320*numel(names) 300]);
    for c = 1:numel(names)
        subplot(1, numel(names), c); hold on;
        for k = 1:numel(methods)
            y = nan(size(values));
            for i = 1:size(rows, 1)
                if strcmp(rows{i, 3}, names{c}) && strcmp(rows{i, 4}, methods{k}) && strcmp(rows{i, 5}, metric)
                    y(values == rows{i, 2}) = rows{i, 6};
                end
            end
            plot(values, y, '-o', 'DisplayName', methods{k});
        end
        title(names{c}, 'Interpreter', 'none'); xlabel(param, 'Interpreter', 'none');
        ylabel(metric, 'Interpreter', 'none'); grid on;
        if c == numel(names), legend('Interpreter', 'none', 'Location', 'best'); end
    end
    kt_save_fig(run, fh, ['fig_sweep_' metric]);
    kt_run_end(run, rows, sprintf('sweep %s over %s (%s plotted)', param, mat2str(values, 4), metric));
catch err
    kt_run_end(run, {}, '', ['FAILED: ' err.message]);
    rethrow(err);
end
end
