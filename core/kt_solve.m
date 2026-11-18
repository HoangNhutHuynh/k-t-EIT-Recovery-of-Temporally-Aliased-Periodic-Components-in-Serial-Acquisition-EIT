function [X, Phi] = kt_solve(rec, Y, tau, nus)
%KT_SOLVE k-t reconstruction of one frequency bin.
%   [X, Phi] = kt_solve(rec, Y, tau, nus)
%   Y   : M x 1 data at the bin (kt_spectrum)
%   tau : M x 1 acquisition offsets (s, kt_timing_model)
%   nus : signed physical frequencies assumed to share this bin;
%         NaN = timing ignored (D = I, conventional reconstruction)
%   Solves
%     min_Z  sum_m w_m |Y_m - sum_l (D_l J Z_l)_m|^2 + lam2 * sum_l Z_l' RtR Z_l,
%     D_l = diag(exp(1i*2*pi*nus(l)*tau)),
%   and returns X_l = Z_l (nus(l) > 0) or conj(Z_l) (nus(l) < 0), i.e. the
%   complex amplitude image of the positive-frequency component.
%   With one frequency the normal matrix equals J'WJ (D is unitary and
%   diagonal), so Z = R*(D'*Y): phase-correct the data, then apply R.
%   Phi is the objective at the optimum (data term in chi-square units).
%
% Author: Hoang Nhut Huynh

J = rec.J;
w = rec.w;
N = size(J, 2);
L = numel(nus);
D = ones(numel(Y), L);
for l = 1:L
    if ~isnan(nus(l))
        D(:, l) = exp(1i*2*pi*nus(l)*tau(:));
    end
end

if L == 1
    Z = rec.R * (conj(D) .* Y(:));
else
    G = zeros(L*N);
    b = zeros(L*N, 1);
    for i = 1:L
        ii = (i-1)*N + (1:N);
        b(ii) = J.' * (w .* conj(D(:, i)) .* Y(:));
        for j = 1:L
            jj = (j-1)*N + (1:N);
            G(ii, jj) = J.' * ((w .* conj(D(:, i)) .* D(:, j)) .* J);
        end
        G(ii, ii) = G(ii, ii) + rec.lam2 * rec.RtR;
    end
    Z = reshape(G \ b, N, L);
end

r = Y(:) - sum((J*Z) .* D, 2);
Phi = real(sum(w .* abs(r).^2) + rec.lam2 * sum(sum(conj(Z) .* (rec.RtR*Z))));
X = Z;
neg = nus < 0;
X(:, neg) = conj(Z(:, neg));
end
