function C = kt_alias_map(nu, fs, tol)
%KT_ALIAS_MAP Where each physical frequency lands after sampling at the frame rate.
%   C = kt_alias_map(nu, fs, tol)
%   nu  : physical frequencies (Hz, > 0)
%   tol : bins closer than tol (Hz) are treated as the same bin
%   For each component c:
%     C(c).nu     physical frequency
%     C(c).bin    analysed (folded) frequency in [0, fs/2]
%     C(c).signed signed physical frequency seen at +bin
%                 (negative = the conjugate image of the component)
%     C(c).ell    alias order: signed = bin + ell*fs
%     C(c).valid  false if the component folds onto DC or Nyquist
%     C(c).group  components sharing a bin (solved jointly by 'kt')
%   Example fs = 44, nu = [10 20 30 40] -> bin [10 20 14 4],
%   signed [10 20 -30 -40], ell [0 0 -1 -1].
%
% Author: Hoang Nhut Huynh

nu = nu(:).';
n = numel(nu);
C = struct('nu', cell(1, n), 'bin', [], 'signed', [], 'ell', [], ...
           'valid', [], 'group', []);
for c = 1:n
    phi = mod(nu(c) + fs/2, fs) - fs/2;   % folded frequency in [-fs/2, fs/2)
    s = 1;
    if phi < 0, s = -1; end
    C(c).nu = nu(c);
    C(c).bin = abs(phi);
    C(c).signed = s*nu(c);
    C(c).ell = round((C(c).signed - C(c).bin)/fs) + 0;   % +0 avoids printing -0
    C(c).valid = C(c).bin > tol && abs(C(c).bin - fs/2) > tol;
end
g = 0;
for c = 1:n
    if ~isempty(C(c).group), continue; end
    g = g + 1;
    for d = c:n
        if isempty(C(d).group) && abs(C(d).bin - C(c).bin) < tol
            C(d).group = g;
        end
    end
end
end
