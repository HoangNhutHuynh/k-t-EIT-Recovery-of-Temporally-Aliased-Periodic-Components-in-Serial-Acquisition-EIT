function q = kt_quadrants(x, cent, area)
%KT_QUADRANTS Fraction of area-weighted |x| in the four image quadrants.
%   q = kt_quadrants(x, cent, area), order:
%   [x<0 & y>=0, x>=0 & y>=0, x<0 & y<0, x>=0 & y<0]
%   (quadrant amplitudes as used by Yerworth et al. 2017).
%
% Author: Hoang Nhut Huynh

v = area(:) .* abs(x(:));
L = cent(:, 1) < 0;
T = cent(:, 2) >= 0;
q = [sum(v(L & T)), sum(v(~L & T)), sum(v(L & ~T)), sum(v(~L & ~T))] / sum(v);
end
