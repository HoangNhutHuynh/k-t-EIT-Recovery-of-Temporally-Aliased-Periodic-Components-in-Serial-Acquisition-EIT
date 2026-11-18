function [mu, lo, hi] = kt_boot_ci(x, B, level)
%KT_BOOT_CI Percentile bootstrap CI of the mean over independent blocks.
%   [mu, lo, hi] = kt_boot_ci(x, B, level)   x: one value per block
%   Set the random seed before calling (rng) for reproducible intervals.
%
% Author: Hoang Nhut Huynh

x = x(:);
x = x(~isnan(x));
n = numel(x);
mu = mean(x);
lo = NaN; hi = NaN;
if n < 2, return; end
bm = sort(mean(x(randi(n, n, B)), 1));
lo = bm(max(1, floor((1 - level)/2*B)));
hi = bm(min(B, ceil((1 + level)/2*B)));
end
