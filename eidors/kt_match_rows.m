function [rows, keyS, keyM] = kt_match_rows(fmdl, pairs)
%KT_MATCH_ROWS Model measurement row of each (injection, measuring pair).
%   rows = kt_match_rows(fmdl, pairs)   pairs(:,1) injection s, pairs(:,2)
%   first electrode m of the adjacent measuring pair (m, m+1).
%   keyS, keyM: injection and first measuring electrode of every model row.
%   Rows follow the fwd_solve output order (injection-major, meas_select applied).
%   The voltage sign convention does not matter: data are used as v/mean(v) - 1.
%
% Author: Hoang Nhut Huynh

nE = size(fmdl.stimulation(1).meas_pattern, 2);
keyS = []; keyM = [];
for s = 1:numel(fmdl.stimulation)
    mp = full(fmdl.stimulation(s).meas_pattern);
    for r = 1:size(mp, 1)
        e = find(mp(r, :) ~= 0);
        if numel(e) ~= 2, error('kt_match_rows: measurement is not a pair'); end
        if mod(e(1), nE) + 1 == e(2), m = e(1); else, m = e(2); end   % first electrode of (m, m+1)
        keyS(end+1, 1) = s; keyM(end+1, 1) = m; %#ok<AGROW>
    end
end
if isfield(fmdl, 'meas_select') && ~isempty(fmdl.meas_select) && numel(fmdl.meas_select) == numel(keyS)
    sel = logical(fmdl.meas_select(:));
    keyS = keyS(sel); keyM = keyM(sel);
end
rows = zeros(size(pairs, 1), 1);
for i = 1:size(pairs, 1)
    j = find(keyS == pairs(i, 1) & keyM == pairs(i, 2));
    if numel(j) ~= 1
        error('kt_match_rows: injection %d, pair (%d,%d) not found in the model', pairs(i,1), pairs(i,2), mod(pairs(i,2), nE) + 1);
    end
    rows(i) = j;
end
end
