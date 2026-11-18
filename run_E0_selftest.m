% run_E0_selftest  Check the k-t core on a synthetic problem (no EIDORS needed).
% Every line must print PASS before running E1-E4.
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;
kt_selftest();
