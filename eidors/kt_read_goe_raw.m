function [v, pairs, info] = kt_read_goe_raw(files, n_elec)
%KT_READ_GOE_RAW Read Goe-MF II .get files directly, keeping device order.
%   [v, pairs, info] = kt_read_goe_raw(files, n_elec)
%   Layout checked on chb_hfov_porcine (CHB_pig06-55.get): 256 float32 per
%   frame; values 1..208 = 16 injections x 13 adjacent measurements, the
%   measurement pairs of injection s starting at electrode s+2 (cyclic);
%   values 209..256 = auxiliary data. Columns that are zero in every frame
%   were not stored (here the reciprocal half: injection s keeps pairs
%   s+2..16, i.e. 104 of 208) and are dropped.
%   v     : M x F voltages of the stored measurements (device order)
%   pairs : M x 3 [injection s, first measuring electrode m, position k]
%   The acquisition time of a row is set by its injection slot s.
%
% Author: Hoang Nhut Huynh

nE = n_elec;
nK = nE - 3;
v = [];
for i = 1:numel(files)
    fid = fopen(files{i}, 'r', 'ieee-le');
    if fid < 0, error('kt_read_goe_raw: cannot open %s', files{i}); end
    x = fread(fid, inf, 'float32=>double');
    fclose(fid);
    if mod(numel(x), 256) ~= 0
        error('kt_read_goe_raw: %s has %d values, not a multiple of 256', files{i}, numel(x));
    end
    x = reshape(x, 256, []);
    v = [v, x(1:nE*nK, :)]; %#ok<AGROW>
    fprintf('read %s (raw float32): %d frames\n', files{i}, size(x, 2));
end
s = repelem((1:nE).', nK);
k = repmat((1:nK).', nE, 1);
m = mod(s + k, nE) + 1;                 % first electrode of the measuring pair
keep = any(v ~= 0, 2);
v = v(keep, :);
pairs = [s(keep), m(keep), k(keep)];
info.format = 'goe_raw';
info.rows_raw = nE*nK;
info.rows_kept = nnz(keep);
info.aux = [];
fprintf('goe_raw: kept %d of %d measurement columns (others never stored)\n', nnz(keep), nE*nK);
end
