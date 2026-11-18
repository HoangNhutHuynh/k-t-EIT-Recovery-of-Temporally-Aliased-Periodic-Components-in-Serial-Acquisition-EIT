function X = kt_try_ftc(dv, fs, nu, tau, rec, opt, dv_full, rows, method)
%KT_TRY_FTC Baseline: EIDORS frame_timing_correct(.., 'FFT'), then timing-free GN.
%   X = kt_try_ftc(dv, fs, nu, tau, rec, opt)
%   X = kt_try_ftc(dv, fs, nu, tau, rec, opt, dv_full, rows)
%   frame_timing_correct (Yerworth and Bayford 2013) assumes 208 rows read one
%   after another over the whole frame. For half-stored data pass the
%   reciprocity-completed frame dv_full (kt_fill_reciprocal, model order) and
%   the model row of every stored measurement (rows); the corrected rows are
%   then taken back to the stored set.
%   Returns [] (and prints why) on any error.
%
% Author: Hoang Nhut Huynh

X = [];
if nargin < 9 || isempty(method), method = 'FFT'; end
if exist('frame_timing_correct', 'file') ~= 2
    fprintf('ftc baseline skipped: frame_timing_correct not found on the path\n');
    return
end
try
    if nargin < 7 || isempty(dv_full)
        dvc = frame_timing_correct(dv, method, []);   % 3 arguments: with 2, EIDORS 3.12 resets method to 'fft'
    else
        dvc = frame_timing_correct(dv_full, method, []);
        dvc = dvc(rows, :);
    end
    if ~isequal(size(dvc), size(dv)), error('output size %s', mat2str(size(dvc))); end
    A = kt_analyze(real(dvc), fs, nu, tau, rec, {'naive'}, opt);
    X = A.X.naive;
catch err
    fprintf('ftc baseline skipped: frame_timing_correct failed: %s\n', err.message);
end
end
