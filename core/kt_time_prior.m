function [X, Xf] = kt_time_prior(rec, dv, fs, bins, win, ts, tw)
%KT_TIME_PRIOR Temporal-prior baseline (Adler, Dai & Lionheart 2007), same J, prior and lambda as k-t.
%   X = kt_time_prior(rec, dv, fs, bins, win, ts, tw)
%   Frame t is reconstructed from frames t-ts..t+ts with a space-time prior
%   whose temporal correlation is tw^|i-j| (data form of EIDORS
%   inv_solve_time_prior), assuming every measurement of a frame is
%   simultaneous. Whitened by W = diag(rec.w), spatial prior P = inv(rec.RtR),
%   noise level rec.lam2:
%     x_t = P Jw' * [ (T kron Jw P Jw' + lam2 I)^-1 * y_stack ]  weighted by row c of T
%   Solved per frame with the eigen-decompositions of T and Jw P Jw'.
%   Component images: windowed DTFT of the frame images at bins.
%   ts = 2, tw = 0.8 as in the EIDORS unit test.
%
% Author: Hoang Nhut Huynh

[M, F] = size(dv);
Jw = sqrt(rec.w) .* rec.J;
P = (rec.RtR + 1e-12*trace(rec.RtR)/size(rec.RtR, 1)*eye(size(rec.RtR))) \ eye(size(rec.RtR));
A = Jw * P * Jw.';
[U, L] = eig((A + A.')/2);
lam = diag(L);
k = (-ts:ts).';
T = tw .^ abs(k - k.');
[V, Th] = eig(T);
th = diag(Th);
c = ts + 1;
Yu = U.' * (sqrt(rec.w) .* dv);                     % M x F
idx = min(max((1:F) + k, 1), F);                    % (2ts+1) x F, edge frames replicated
Q = reshape(Yu(:, idx(:)), M, numel(k), F);         % M x l x F
Qm = reshape(permute(Q, [1 3 2]), M*F, numel(k));   % rows: m + (f-1)M
Qm = Qm * V;
Qm = Qm ./ repmat(lam * th.' + rec.lam2, F, 1);
Wt = reshape(Qm * (V.' * T(c, :).'), M, F);
Xf = P * (Jw.' * (U * Wt));                          % N x F frame images
X = kt_spectrum(Xf, fs, bins, win);
end
