function A = kt_analyze(dv, fs, nu, tau, rec, methods, opt)
%KT_ANALYZE Component images for every method from one multi-frame record.
%   A = kt_analyze(dv, fs, nu, tau, rec, methods, opt)
%   dv  : M x F normalised difference data
%   nu  : physical frequencies of the components to image (Hz)
%   tau : M x 1 acquisition offsets assumed by the timing-aware methods
%   rec : kt_recon_matrix output (same R, lam2 for all methods)
%   methods (cell of char):
%     'naive'     timing ignored (= per-frame GN, then projection on the bin)
%     'yerworth'  phase correction at the folded frequency (ell = 0),
%                 i.e. the FFT timing correction of Yerworth et al. 2017
%     'kt_single' phase correction at the true alias frequency, one component
%     'kt'        as kt_single, but components sharing a bin are solved jointly
%     'kt_auto'   alias order chosen by the smallest MAP objective,
%                 |ell| <= opt.l_max (no knowledge of the harmonic structure)
%   opt.window ('hann'|'rect'), opt.tol (Hz), opt.l_max
%   Output: A.C alias map, A.Y (M x C), A.X.(method) (N x C complex),
%           A.Phi.(method) (1 x C), A.nu_used.(method) (signed frequency used)
%
% Author: Hoang Nhut Huynh

C = kt_alias_map(nu, fs, opt.tol);
nC = numel(C);
Y = kt_spectrum(dv, fs, [C.bin], opt.window);
N = size(rec.J, 2);
A.C = C;
A.Y = Y;
for k = 1:numel(methods)
    m = methods{k};
    X = nan(N, nC);
    Phi = nan(1, nC);
    used = nan(1, nC);
    for c = 1:nC
        if ~C(c).valid, continue; end
        y = Y(:, c);
        switch m
            case 'naive'
                [x, p] = kt_solve(rec, y, tau, NaN);
            case 'yerworth'
                [x, p] = kt_solve(rec, y, tau, C(c).bin);
                used(c) = C(c).bin;
            case 'kt_single'
                [x, p] = kt_solve(rec, y, tau, C(c).signed);
                used(c) = C(c).signed;
            case 'kt'
                grp = find([C.group] == C(c).group & [C.valid]);
                [xg, p] = kt_solve(rec, y, tau, [C(grp).signed]);
                x = xg(:, grp == c);
                used(c) = C(c).signed;
            case 'kt_auto'
                p = inf;
                for ell = -opt.l_max:opt.l_max
                    nu_try = C(c).bin + ell*fs;
                    if abs(nu_try) < opt.tol, continue; end
                    if ell ~= 0 && isfield(opt, 'nu_max') && numel(opt.nu_max) >= c && abs(nu_try) > opt.nu_max(c), continue; end
                    [xt, pt] = kt_solve(rec, y, tau, nu_try);
                    if pt < p
                        p = pt; x = xt; used(c) = nu_try;
                    end
                end
            otherwise
                error('kt_analyze: unknown method "%s"', m);
        end
        X(:, c) = x;
        Phi(c) = p;
    end
    A.X.(m) = X;
    A.Phi.(m) = Phi;
    A.nu_used.(m) = used;
end
end
