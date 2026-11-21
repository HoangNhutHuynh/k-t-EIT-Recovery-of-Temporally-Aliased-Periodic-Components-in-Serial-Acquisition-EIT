function kt_write_rows(fid, tag, rec, rows)
%KT_WRITE_ROWS Append result rows {param, pvalue, component, method, metric, value} to an open CSV.
%   Columns written: tag, record, param, pvalue, component, method, metric, value
%
% Author: Hoang Nhut Huynh

for i = 1:size(rows, 1)
    fprintf(fid, '%s,%s,%s,%g,%s,%s,%s,%.6g\n', tag, rec, char(string(rows{i, 1})), rows{i, 2}, ...
            char(string(rows{i, 3})), char(string(rows{i, 4})), char(string(rows{i, 5})), rows{i, 6});
end
end
