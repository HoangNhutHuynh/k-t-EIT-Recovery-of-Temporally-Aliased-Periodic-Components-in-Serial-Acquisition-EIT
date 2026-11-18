function out = kt_separability(J, w, tau, dnu, ranks)
%KT_SEPARABILITY Principal angles between the data subspaces of two aliases.
%   out = kt_separability(J, w, tau, dnu, ranks)
%   Two components at signed frequencies nu_a and nu_b = nu_a + dnu share a
%   bin. With W^(1/2) J = U S V', the measurable data subspaces are
%   span(D_a U_r) and span(D_b U_r). The cosines of their principal angles are
%   the singular values of  U_r' * diag(exp(1i*2*pi*dnu*tau)) * U_r.
%   Angle 0 = inseparable (always the case for simultaneous sampling, tau = 0).
%   Larger rank r = finer spatial detail kept = smaller angles: this is the
%   space-time trade-off of k-t EIT.
%   out.min_angle_deg(k), out.mean_cos2(k) for r = ranks(k).
%
% Author: Hoang Nhut Huynh

[U, ~, ~] = svd(sqrt(w(:)) .* J, 'econ');
ph = exp(1i*2*pi*dnu*tau(:));
ranks = ranks(ranks <= size(U, 2));
out.ranks = ranks;
out.min_angle_deg = zeros(size(ranks));
out.mean_cos2 = zeros(size(ranks));
for k = 1:numel(ranks)
    Ur = U(:, 1:ranks(k));
    s = min(svd(Ur' * (ph .* Ur)), 1);
    out.min_angle_deg(k) = acos(max(s)) * 180/pi;
    out.mean_cos2(k) = mean(s.^2);
end
end
