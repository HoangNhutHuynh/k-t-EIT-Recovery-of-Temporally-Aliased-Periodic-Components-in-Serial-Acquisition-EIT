function lines = kt_flatten(s, prefix)
%KT_FLATTEN Every field of a (nested) struct as 'a.b.c = value' lines.
%
% Author: Hoang Nhut Huynh

lines = {};
f = fieldnames(s);
for i = 1:numel(f)
    v = s.(f{i});
    name = [prefix f{i}];
    if isstruct(v) && isscalar(v)
        lines = [lines, kt_flatten(v, [name '.'])]; %#ok<AGROW>
    elseif isstruct(v)
        for j = 1:numel(v)
            lines = [lines, kt_flatten(v(j), sprintf('%s(%d).', name, j))]; %#ok<AGROW>
        end
    else
        lines{end+1} = sprintf('%s = %s', name, kt_val2str(v)); %#ok<AGROW>
    end
end
end
