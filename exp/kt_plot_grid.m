function fh = kt_plot_grid(I, imgs, colnames, rownames, own_scale, ttl, ref_col, opts)
%KT_PLOT_GRID Component images (rows) for every method (columns), IEEE two-column layout.
%   fh = kt_plot_grid(I, imgs, colnames, rownames, own_scale, ttl, ref_col, opts)
%   imgs      : cell, one N x C complex matrix per column
%   own_scale : logical per column; true = column scaled to its own maximum
%               (e.g. phantom truth); false = scaled by the row maximum of the
%               reference column ref_col (or of all such columns); values beyond
%               the scale saturate
%   ttl       : stored as the figure name only (put titles in the caption)
%   opts.mode : 'abs'     |X|, sequential colour map (default)
%               'aligned' Re(X e^{i theta}) after the global phase (and frequency
%                         sign) that best matches the reference column, diverging
%                         colour map. This is what rho_ref measures, so phase
%                         errors that |X| hides become visible.
%               'error'   reference column: |X_ref|; other columns: |X e^{i theta} - X_ref|
%                         (same alignment as 'aligned'), white-to-red map, so where a
%                         method departs from the reference (e.g. inside the lungs)
%                         is visible directly; panels print e (relative error) and
%                         L (the same inside the anatomical lung mask, if available)
%   opts.rho  : true -> print rho against the reference column in each panel
%   Display: element values averaged to the nodes (area weighted), interpolated
%   shading; anatomical lung outline (I.lung_poly) and electrode 1 (sternum,
%   triangle) overlaid; one shared colour bar. Width 7.16 in, 8 pt text.
%
% Author: Hoang Nhut Huynh

if nargin < 7, ref_col = []; end
if nargin < 8, opts = struct(); end
if ~isfield(opts, 'mode'), opts.mode = 'abs'; end
if ~isfield(opts, 'rho'), opts.rho = ~isempty(ref_col); end
aligned = any(strcmp(opts.mode, {'aligned', 'error'})) && ~isempty(ref_col);
errmode = strcmp(opts.mode, 'error') && ~isempty(ref_col);
nR = numel(rownames);
nK = numel(imgs);
W = 7.16;                                  % in
lm = 0.30; rm = 0.62; tm = 0.22; bm = 0.05;
cw = (W - lm - rm) / nK;
ch = cw;
H = tm + bm + nR*ch;
fh = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 W H], 'Name', ttl, ...
            'PaperUnits', 'inches', 'PaperPosition', [0 0 W H], 'PaperSize', [W H]);
try theme(fh, 'light'); catch, end
if errmode
    t = linspace(0, 1, 256).';
    cmap = interp1([0 0.35 0.7 1], [1 1 1; 0.99 0.75 0.45; 0.89 0.29 0.20; 0.55 0.0 0.05], t);
    clim = [0 1]; olc = [0 0 0]; cbl = '|X_{ref}|  and  |X - X_{ref}|, / max|X_{ref}|';
elseif aligned
    t = linspace(0, 1, 128).';
    cmap = [[0.23 + 0.77*t, 0.30 + 0.70*t, 0.75 + 0.25*t]; ...
            [1 - 0.29*t, 1 - 0.98*t, 1 - 0.85*t]];          % blue - white - red
    clim = [-1 1]; olc = [0 0 0]; cbl = 'Re(X) / max|ref|';
else
    try cmap = parula(256); catch, cmap = viridis(256); end
    clim = [0 1]; olc = [1 1 1]; cbl = '|X| / max|ref|';
end
set(fh, 'Colormap', cmap);

Ne = size(I.elems, 1); Nn = size(I.nodes, 1);
Aen = sparse(repmat((1:Ne).', 3, 1), I.elems(:), repmat(I.area(:), 3, 1), Ne, Nn);
den = full(sum(Aen, 1)).';
den(den == 0) = 1;
nodal = @(v) (Aen.' * v(:)) ./ den;
a = I.area(:);
p1 = mean(I.nodes(I.fmdl.electrode(1).nodes, 1:2), 1);
E = sort([I.elems(:, [1 2]); I.elems(:, [2 3]); I.elems(:, [3 1])], 2);   % boundary edges (used once)
[Eu, ~, ie] = unique(E, 'rows');
Eb = Eu(accumarray(ie, 1) == 1, :);
bx = [I.nodes(Eb(:, 1), 1), I.nodes(Eb(:, 2), 1), nan(size(Eb, 1), 1)].';
by = [I.nodes(Eb(:, 1), 2), I.nodes(Eb(:, 2), 2), nan(size(Eb, 1), 1)].';
R = max(sqrt(sum(I.nodes.^2, 2)));
lim = 1.08 * R * [-1 1];

for r = 1:nR
    if aligned
        z = imgs{ref_col}(:, r);
        z = z * exp(-0.5i * angle(sum(a .* z.^2)));          % reference: most real global phase
    end
    vals = cell(1, nK); ex = cell(1, nK);
    for k = 1:nK
        x = imgs{k}(:, r);
        if aligned
            if k == ref_col
                x = z;
            else
                p = sum(a .* conj(x) .* z); q = sum(a .* x .* z);
                if abs(q) > abs(p), x = conj(x); p = sum(a .* conj(x) .* z); end
                x = x * exp(1i * angle(p));
            end
            if errmode
                if k == ref_col, vals{k} = nodal(abs(z)); else, vals{k} = nodal(abs(x - z)); end
                ex{k} = x - z;
            else
                vals{k} = nodal(real(x));
            end
        else
            vals{k} = nodal(abs(x));
        end
    end
    if ~isempty(ref_col)
        shared = max(abs(vals{ref_col}));
    else
        shared = 0;
        for k = 1:nK, if ~own_scale(k), shared = max(shared, max(abs(vals{k}))); end, end
    end
    for k = 1:nK
        top = shared;
        if own_scale(k) || ~(top > 0), top = max(abs(vals{k})); end
        if ~(top > 0), top = 1; end
        vn = max(min(vals{k} / top, clim(2)), clim(1));
        ax = axes('Parent', fh, 'Units', 'inches', ...
                  'Position', [lm + (k-1)*cw, H - tm - r*ch, cw*0.98, ch*0.98]);
        patch(ax, 'Faces', I.elems, 'Vertices', I.nodes, 'FaceVertexCData', vn, ...
              'FaceColor', 'interp', 'EdgeColor', 'none');
        hold(ax, 'on');
        if isfield(I, 'lung_poly')
            for q = 1:numel(I.lung_poly)
                L = I.lung_poly{q};
                plot(ax, L([1:end 1], 1), L([1:end 1], 2), '-', 'Color', olc, 'LineWidth', 0.5 + 0.4*errmode);
            end
        end
        plot(ax, bx(:), by(:), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
        plot(ax, p1(1)*1.04, p1(2)*1.04, 'v', 'MarkerSize', 3, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
        axis(ax, 'equal'); set(ax, 'XLim', lim, 'YLim', lim, 'Visible', 'off');
        caxis(ax, clim);
        if r == 1
            text(ax, 0, lim(2), colnames{k}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                 'FontSize', 8, 'FontWeight', 'bold', 'Interpreter', 'none', 'Color', 'k');
        end
        if k == 1
            text(ax, lim(1), 0, rownames{r}, 'Rotation', 90, 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'bottom', 'FontSize', 8, 'Interpreter', 'none', 'Color', 'k');
        end
        if errmode && k ~= ref_col
            e = sqrt(sum(a .* abs(ex{k}).^2) / sum(a .* abs(z).^2));
            txt = sprintf('e=%.2f', e);
            if isfield(I, 'lung')
                l = I.lung(:);
                txt = sprintf('%s, L=%.2f', txt, sqrt(sum(a(l) .* abs(ex{k}(l)).^2) / sum(a(l) .* abs(z(l)).^2)));
            end
            text(ax, lim(2), lim(1), txt, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                 'FontSize', 7, 'Color', 'k');
        elseif opts.rho && ~isempty(ref_col) && k ~= ref_col
            s = kt_compare(imgs{k}(:, r), imgs{ref_col}(:, r), a);
            text(ax, lim(2), lim(1), sprintf('\\rho=%.2f', s.rho), 'HorizontalAlignment', 'right', ...
                 'VerticalAlignment', 'bottom', 'FontSize', 7, 'Color', 'k');
        end
    end
end
cpos = [W - rm + 0.08, bm + 0.1, 0.1, nR*ch - 0.2];
cax = axes('Parent', fh, 'Units', 'inches', 'Position', cpos, 'Visible', 'off');
set(cax, 'Colormap', cmap);
caxis(cax, clim);
cb = colorbar(cax);
set(cb, 'Units', 'inches', 'Position', cpos, 'FontSize', 7, 'Color', [0.15 0.15 0.15]);
try set(cb, 'Ticks', unique([clim(1), mean(clim), clim(2)])); catch, end
ylabel(cb, cbl, 'FontSize', 7, 'Color', 'k');
end
