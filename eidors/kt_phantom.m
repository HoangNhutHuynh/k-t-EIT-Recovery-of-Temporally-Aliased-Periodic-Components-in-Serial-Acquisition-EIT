function [nu, X, names, analyze] = kt_phantom(cfg, cent, R)
%KT_PHANTOM Complex amplitude images of the simulated components.
%   [nu, X, names, analyze] = kt_phantom(cfg, cent, R)
%   Conductivity change of element e at time t:
%     delta_e(t) = sum_c Re( X(e,c) * exp(1i*2*pi*nu(c)*t) )
%   Regions are analytic, so the truth is evaluated on any mesh (fine forward
%   mesh for the data, coarse inverse mesh for scoring).
%
% Author: Hoang Nhut Huynh

comp = cfg.sim.components;
nC = numel(comp);
nu = zeros(1, nC);
X = zeros(size(cent, 1), nC);
names = {comp.name};
analyze = [comp.analyze];
for c = 1:nC
    nu(c) = comp(c).nu;
    if comp(c).harm > 0
        nu(c) = comp(c).harm * cfg.sim.f0;
    end
    lr = exp(1i*comp(c).phase_lr_deg*pi/180 * (cent(:, 1) > 0));   % phase lag of the x > 0 half
    X(:, c) = comp(c).amp * exp(1i*comp(c).phase_deg*pi/180) * kt_region(comp(c).region, cent, R) .* lr;
end
end

function mask = kt_region(name, cent, R)
% Circular cross-section of radius R, image orientation: x < 0 = patient right.
x = cent(:, 1)/R; y = cent(:, 2)/R;
rlung = ((x + 0.45)/0.28).^2 + ((y - 0.05)/0.45).^2 <= 1;
llung = ((x - 0.45)/0.28).^2 + ((y - 0.05)/0.45).^2 <= 1;
heart = (x/0.15).^2 + ((y + 0.05)/0.15).^2 <= 1;
switch name
    case 'lungs',      mask = rlung | llung;
    case 'right_lung', mask = rlung;
    case 'left_lung',  mask = llung;
    case 'heart',      mask = heart;
    case 'global',     mask = true(size(x));
    otherwise, error('kt_phantom: unknown region "%s"', name);
end
mask = double(mask);
end
