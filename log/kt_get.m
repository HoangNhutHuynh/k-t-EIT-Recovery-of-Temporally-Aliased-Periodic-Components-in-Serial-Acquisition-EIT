function v = kt_get(cfg, path, default)
%KT_GET Value of a dotted config path, e.g. kt_get(cfg, 'rec.alpha').
%   v = kt_get(cfg, path, default) returns default when the path does not exist
%   (configs saved before a field was added).
%
% Author: Hoang Nhut Huynh

parts = strsplit(path, '.');
if nargin < 3
    v = getfield(cfg, parts{:});
    return
end
v = default;
s = cfg;
for k = 1:numel(parts)
    if ~isstruct(s) || ~isfield(s, parts{k}), return; end
    s = s.(parts{k});
end
v = s;
end
