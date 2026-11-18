function [v, info] = kt_load_real(Rcfg, n_elec)
%KT_LOAD_REAL Read raw EIT voltages (M x F) with eidors_readdata.
%   [v, info] = kt_load_real(cfg.real, cfg.acq.n_elec)
%   Tries the formats in Rcfg.formats in order and logs which one worked.
%   Measurement order matters for k-t EIT: 'GET-RAW' should keep the device
%   order, 'GET' may reorder (see eidors_readdata help). A remaining cyclic
%   offset / reversal of the injection order is estimated by kt_estimate_timing.
%   256 rows (all adjacent pairs) are reduced to the 208 that do not touch
%   the current electrodes, assuming injection-major order.
%
% Author: Hoang Nhut Huynh

v = [];
info.files = Rcfg.files;
info.format = '';
info.rows_raw = NaN;
info.aux = [];
for i = 1:numel(Rcfg.files)
    fname = Rcfg.files{i};
    if ~exist(fname, 'file')
        error('kt_load_real: file not found: %s (edit cfg.real.files)', fname);
    end
    vi = [];
    for k = 1:numel(Rcfg.formats)
        try
            [vi, aux] = eidors_readdata(fname, Rcfg.formats{k});
            info.format = Rcfg.formats{k};
            info.aux = aux;
            break
        catch err
            fprintf('eidors_readdata(%s, %s) failed: %s\n', fname, Rcfg.formats{k}, err.message);
        end
    end
    if isempty(vi), error('kt_load_real: no format could read %s', fname); end
    info.rows_raw = size(vi, 1);
    fprintf('read %s as %s: %d rows x %d frames\n', fname, info.format, size(vi, 1), size(vi, 2));
    v = [v, vi]; %#ok<AGROW>
end
if isstruct(info.aux) && isfield(info.aux, 'frame_rate')
    fprintf('header frame_rate = %g Hz (cfg.real.fs = %g Hz is used)\n', info.aux.frame_rate, Rcfg.fs);
end
nE = n_elec;
if size(v, 1) == nE*nE
    keep = true(nE, nE);              % (meas pair, injection pair), injection-major
    for s = 1:nE
        sp = [s, mod(s, nE) + 1];
        for m = 1:nE
            mp = [m, mod(m, nE) + 1];
            keep(m, s) = ~any(ismember(mp, sp));
        end
    end
    v = v(keep(:), :);
    fprintf('reduced %d rows to %d (no measurement on current electrodes)\n', nE*nE, size(v, 1));
elseif size(v, 1) ~= nE*(nE - 3)
    error('kt_load_real: %d rows per frame, expected %d or %d', size(v, 1), nE*nE, nE*(nE - 3));
end
v = real(v);
end
