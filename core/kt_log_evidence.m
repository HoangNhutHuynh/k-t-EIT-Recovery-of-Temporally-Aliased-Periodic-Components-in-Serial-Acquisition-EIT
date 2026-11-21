function nle = kt_log_evidence(rec, y, tau, nus, w)
%KT_LOG_EVIDENCE Negative log marginal likelihood of one frequency bin (R3 review).
%   nle = kt_log_evidence(rec, y, tau, nus)
%   nle = kt_log_evidence(rec, y, tau, nus, w)   noise weights w override rec.w
%   The objective of kt_solve,  sum_m w_m |r_m|^2 + lam2 * sum_l Z_l' RtR Z_l,
%   is the negative log posterior of the complex Gaussian model
%       y = sum_l D_l J Z_l + n,   n ~ CN(0, diag(1./w)),   Z_l ~ CN(0, inv(lam2*RtR)),
%   so y ~ CN(0, C) with  C = diag(1./w) + sum_l D_l J inv(lam2*RtR) J' D_l'.
%   nle = y' inv(C) y + log det C   (the constant M*log(pi) is dropped).
%   y' inv(C) y equals the minimised objective; log det C charges the extra unknowns of a
%   larger model, so nle can compare models with different numbers of components (for
%   example with and without an unaliased line at the folded bin of an aliased component)
%   without knowing which bins collide. Lower nle = more probable model.
%   nus: signed frequencies as in kt_solve (NaN = timing ignored).
%
% Author: Hoang Nhut Huynh

if nargin < 5 || isempty(w), w = rec.w; end
J = rec.J;
M = size(J, 1);
P = rec.lam2 * rec.RtR;
P = P + 1e-10 * trace(P) / size(P, 1) * eye(size(P));   % guards a singular prior (Laplace)
Q = J * (P \ J.');                                       % M x M prior covariance of J*Z
C = diag(1 ./ w(:));
for l = 1:numel(nus)
    if isnan(nus(l)), d = ones(M, 1); else, d = exp(1i*2*pi*nus(l)*tau(:)); end
    C = C + (d .* Q) .* d';                              % D Q D^H
end
C = (C + C') / 2;
L = chol(C, 'lower');
q = L \ y(:);
nle = real(q' * q) + 2 * sum(log(real(diag(L))));
end
