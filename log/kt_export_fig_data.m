function kt_export_fig_data(figfile)
%KT_EXPORT_FIG_DATA Write the plotted data of a kt_plot_grid figure (.fig) to a plain .mat.
%   kt_export_fig_data('runs/<run>/fig_decim_error.fig') -> fig_decim_error_data.mat
%   For every image panel: axes position (inches), mesh (Vertices, Faces), the plotted
%   nodal values (FaceVertexCData, already scaled to the colour limits), outline lines
%   and panel texts. Used to redraw the manuscript figures outside MATLAB with a
%   consistent IEEE style; no reconstruction is repeated.
%
% Author: Hoang Nhut Huynh

fh = openfig(figfile, 'invisible');
cleanup = onCleanup(@() close(fh));
axs = findobj(fh, 'Type', 'axes');
P = struct('pos', {}, 'V', {}, 'F', {}, 'C', {}, 'clim', {}, 'lx', {}, 'ly', {}, 'lc', {}, 'lw', {}, ...
           'mx', {}, 'my', {}, 'txt', {}, 'txtpos', {});
for i = 1:numel(axs)
    ax = axs(i);
    pa = findobj(ax, 'Type', 'patch');
    if isempty(pa), continue; end
    set(ax, 'Units', 'inches');
    p.pos = get(ax, 'Position');
    p.V = pa(1).Vertices; p.F = pa(1).Faces; p.C = pa(1).FaceVertexCData;
    p.clim = get(ax, 'CLim');
    ln = findobj(ax, 'Type', 'line');
    p.lx = {}; p.ly = {}; p.lc = []; p.lw = []; p.mx = []; p.my = [];
    for j = 1:numel(ln)
        if ~strcmp(ln(j).LineStyle, 'none')
            p.lx{end+1} = ln(j).XData; p.ly{end+1} = ln(j).YData;
            p.lc(end+1, :) = ln(j).Color; p.lw(end+1) = ln(j).LineWidth;
        else
            p.mx = ln(j).XData; p.my = ln(j).YData;      % electrode-1 marker
        end
    end
    tx = findobj(ax, 'Type', 'text');
    p.txt = {}; p.txtpos = zeros(0, 3);
    for j = 1:numel(tx)
        s = tx(j).String; if iscell(s), s = strjoin(s, ' '); end
        if isempty(s), continue; end
        p.txt{end+1} = s; p.txtpos(end+1, :) = [tx(j).Position(1:2), tx(j).Rotation];
    end
    P(end+1) = p; %#ok<AGROW>
end
cm = get(fh, 'Colormap');
[d, n] = fileparts(figfile);
save(fullfile(d, [n '_data.mat']), 'P', 'cm', '-v7');
fprintf('%s: %d panels -> %s\n', n, numel(P), fullfile(d, [n '_data.mat']));
end
