function kt_light_fig(fh)
%KT_LIGHT_FIG Force a white background and dark text on a figure before export.
%   MATLAB R2025a+ follows the desktop theme: in dark mode figures, axes,
%   legends and text are drawn light-on-black and saved that way. This sets
%   the figure theme to light (when available) and then fixes any remaining
%   dark backgrounds or light text explicitly.
%
% Author: Hoang Nhut Huynh

try
    theme(fh, 'light');
catch
end
dark = [0.15 0.15 0.15];
set(fh, 'Color', 'w');
try set(fh, 'InvertHardcopy', 'off'); catch, end
ax = findall(fh, 'Type', 'axes');
for a = ax(:).'
    if strcmp(get(a, 'Visible'), 'on')
        set(a, 'Color', 'w', 'XColor', dark, 'YColor', dark, 'ZColor', dark, 'GridColor', dark);
    end
    for h = [get(a, 'Title'), get(a, 'XLabel'), get(a, 'YLabel'), get(a, 'ZLabel')]
        set(h, 'Color', 'k');
    end
end
lg = findall(fh, 'Type', 'legend');
set(lg, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.6 0.6 0.6]);
cb = findall(fh, 'Type', 'colorbar');
set(cb, 'Color', dark);
tx = [findall(fh, 'Type', 'text'); findall(fh, 'Type', 'textboxshape')];
for t = tx(:).'
    c = get(t, 'Color');
    if isnumeric(c) && numel(c) == 3 && sum(c) > 2.4      % light text (dark theme) -> black
        set(t, 'Color', 'k');
    end
end
end
