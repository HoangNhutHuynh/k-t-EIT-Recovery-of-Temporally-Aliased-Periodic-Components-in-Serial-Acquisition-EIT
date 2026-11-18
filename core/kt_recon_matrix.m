function rec = kt_recon_matrix(J, w, prior, alpha, RtR_ext)
%KT_RECON_MATRIX One-step Gauss-Newton (Tikhonov) matrix shared by all methods.
%   rec = kt_recon_matrix(J, w, prior, alpha, RtR_ext)
%   R    = (J'WJ + lam2*RtR) \ J'W,   W = diag(w)
%   lam2 = alpha * trace(J'WJ) / trace(RtR)   (alpha is dimensionless)
%   prior: 'noser'    RtR = diag(J'WJ)
%          'tikhonov' RtR = I
%          'external' RtR = RtR_ext (e.g. Laplace from EIDORS calc_RtR_prior)
%   Every method uses the same R and lam2, so differences between methods
%   come only from the data model (timing / alias frequency).
%
% Author: Hoang Nhut Huynh

if nargin < 5, RtR_ext = []; end
w = w(:);
JW = J.' .* w.';                 % N x M
H = JW * J;
switch prior
    case 'noser',    RtR = diag(diag(H));
    case 'tikhonov', RtR = eye(size(H));
    case 'external', RtR = full(RtR_ext);
    otherwise, error('kt_recon_matrix: unknown prior "%s"', prior);
end
lam2 = alpha * trace(H) / trace(RtR);
rec.J = J;
rec.w = w;
rec.H = H;
rec.RtR = RtR;
rec.lam2 = lam2;
rec.R = (H + lam2*RtR) \ JW;
end
