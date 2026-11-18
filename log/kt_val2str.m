function s = kt_val2str(v)
%KT_VAL2STR Short text form of a config value (for logs).
%
% Author: Hoang Nhut Huynh

if ischar(v)
    s = ['''' v ''''];
elseif isnumeric(v) || islogical(v)
    if numel(v) <= 40
        s = mat2str(v, 6);
    else
        s = sprintf('[%s %s, %g..%g]', mat2str(size(v)), class(v), min(v(:)), max(v(:)));
    end
elseif iscell(v)
    parts = cellfun(@kt_val2str, v, 'UniformOutput', false);
    s = ['{' strjoin(parts(:).', ' ') '}'];
elseif isa(v, 'function_handle')
    s = func2str(v);
elseif isstruct(v)
    s = sprintf('<struct %s>', mat2str(size(v)));
else
    s = sprintf('<%s>', class(v));
end
end
