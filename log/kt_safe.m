function s = kt_safe(s)
%KT_SAFE Text safe for file names and unquoted CSV fields.
%
% Author: Hoang Nhut Huynh

if isempty(s), s = '-'; return; end
s = regexprep(char(s), '[^A-Za-z0-9_.+\-]', '-');
end
