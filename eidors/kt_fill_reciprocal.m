function dvf = kt_fill_reciprocal(dv, pairs, fmdl)
%KT_FILL_RECIPROCAL Complete a half-stored adjacent data set by reciprocity.
%   dvf = kt_fill_reciprocal(dv, pairs, fmdl)
%   dv    : M x F stored data (device order), pairs from kt_read_goe_raw
%   dvf   : all model rows (fwd_solve order). A missing row (inject s, measure m)
%           is copied from the stored row (inject m, measure s): transfer
%           impedances are reciprocal, and v./mean(v) - 1 removes any sign.
%   Used only to run baselines that need the full 208-row frame (EIDORS
%   frame_timing_correct); the copies carry the acquisition time of the
%   stored measurement, which such baselines ignore anyway.
%
% Author: Hoang Nhut Huynh

[~, keyS, keyM] = kt_match_rows(fmdl, pairs(1, :));
dvf = nan(numel(keyS), size(dv, 2));
for r = 1:numel(keyS)
    i = find(pairs(:, 1) == keyS(r) & pairs(:, 2) == keyM(r), 1);
    if isempty(i), i = find(pairs(:, 1) == keyM(r) & pairs(:, 2) == keyS(r), 1); end
    if isempty(i)
        error('kt_fill_reciprocal: no stored or reciprocal measurement for injection %d, pair %d', keyS(r), keyM(r));
    end
    dvf(r, :) = dv(i, :);
end
end
