function kt_addpaths()
%KT_ADDPATHS Add the k-t EIT folders to the MATLAB path.
%
% Author: Hoang Nhut Huynh

r = fileparts(mfilename('fullpath'));
addpath(r, fullfile(r, 'core'), fullfile(r, 'eidors'), fullfile(r, 'log'), fullfile(r, 'exp'));
end
