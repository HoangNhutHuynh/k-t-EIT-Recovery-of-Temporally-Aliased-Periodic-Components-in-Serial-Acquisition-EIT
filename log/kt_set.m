function cfg = kt_set(cfg, path, v)
%KT_SET Set a dotted config path, e.g. cfg = kt_set(cfg, 'rec.alpha', 0.1).
%   The path must already exist in kt_config (protects against typos).
%
% Author: Hoang Nhut Huynh

parts = strsplit(path, '.');
kt_get(cfg, path);               % errors if the field does not exist
cfg = setfield(cfg, parts{:}, v);
end
