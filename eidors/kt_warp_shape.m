function [fmdl, lungs, info] = kt_warp_shape(fmdl, name, rot_deg)
%KT_WARP_SHAPE Map a circular 2-D EIDORS model onto a thorax outline from shape_library.
%   [fmdl, lungs, info] = kt_warp_shape(fmdl, name, rot_deg)
%   Every node at polar (r, theta) moves radially to r * rho(theta), where
%   rho(theta) is the radius of the normalised library boundary; the mesh
%   topology and the electrode node lists are kept, so the measurement order
%   is unchanged. The outline is rotated (and mirrored if the numbering
%   directions differ) so that library electrode 1 lands on model electrode 1.
%   Libraries without electrode positions (e.g. 'neonate') are assumed to have
%   electrode 1 at the top (+y, sternum) numbered clockwise, as for Goe-MF II
%   data. rot_deg adds a further rotation (degrees).
%   lungs : cell of [x y] polygons in model coordinates (fields whose name
%           contains 'lung', excluding 'lungs_and_heart').
%
% Author: Hoang Nhut Huynh

S = shape_library('get', name);
B = S.boundary(:, 1:2);
c0 = mean(B, 1);
B = B - c0;
sc = max(sqrt(sum(B.^2, 2)));
B = B / sc;
f = fieldnames(S);
lungs = {};
for k = 1:numel(f)
    if ~isempty(strfind(lower(f{k}), 'lung')) && isempty(strfind(lower(f{k}), 'heart')) && isnumeric(S.(f{k}))
        lungs{end+1} = (S.(f{k})(:, 1:2) - c0) / sc; %#ok<AGROW>
    end
end
if isfield(S, 'electrodes') && size(S.electrodes, 1) >= 2
    E = (S.electrodes(:, 1:2) - c0) / sc;
    thL = atan2(E(1:2, 2), E(1:2, 1));
else
    thL = [pi/2; pi/2 - 2*pi/16];                % electrode 1 at the top, clockwise
end
thM = zeros(2, 1);
for e = 1:2
    p = mean(fmdl.nodes(fmdl.electrode(e).nodes, 1:2), 1);
    thM(e) = atan2(p(2), p(1));
end
dirL = sign(wrap(thL(2) - thL(1)));
dirM = sign(wrap(thM(2) - thM(1)));
mirror = dirL ~= dirM;
if mirror
    B(:, 1) = -B(:, 1);
    lungs = cellfun(@(L) [-L(:, 1), L(:, 2)], lungs, 'UniformOutput', false);
    thL = atan2(sin(thL), -cos(thL));
end
phi = thM(1) - thL(1) + rot_deg*pi/180;
Rm = [cos(phi), -sin(phi); sin(phi), cos(phi)];
B = (Rm * B.').';
lungs = cellfun(@(L) (Rm * L.').', lungs, 'UniformOutput', false);

[tb, rb] = cart2pol(B(:, 1), B(:, 2));
tp = [tb - 2*pi; tb; tb + 2*pi];
rp = [rb; rb; rb];
[tp, iu] = unique(tp);
rp = rp(iu);
R0 = max(sqrt(sum(fmdl.nodes(:, 1:2).^2, 2)));
[th, r] = cart2pol(fmdl.nodes(:, 1), fmdl.nodes(:, 2));
rho = interp1(tp, rp, th, 'linear');
fmdl.nodes(:, 1:2) = [r .* rho .* cos(th), r .* rho .* sin(th)];
lungs = cellfun(@(L) L * R0, lungs, 'UniformOutput', false);
info.mirror = mirror;
info.rotation_deg = phi*180/pi;
info.electrode1_deg = thM(1)*180/pi;
fprintf('shape %s: rotated %.1f deg, mirrored %d, %d lung outline(s)\n', name, info.rotation_deg, mirror, numel(lungs));
end

function a = wrap(a)
a = mod(a + pi, 2*pi) - pi;
end
