function m = kt_compare(xh, x, area)
%KT_COMPARE Similarity of a complex component image to a reference image.
%   m = kt_compare(xh, x, area)
%   Invariant to a global phase and to complex conjugation (the sign of the
%   frequency), so baselines are not penalised for an arbitrary convention.
%   m.rho    = |<xh,x>| / (||xh|| ||x||)        1 = same complex pattern
%   m.amp    = ||xh|| / ||x||                   1 = same amplitude
%   m.relerr = min_theta ||xh e^{i theta} - x|| / ||x||
%   Inner products are area weighted (area = element areas).
%
% Author: Hoang Nhut Huynh

if nargin < 3, area = ones(size(x)); end
a = area(:); xh = xh(:); x = x(:);
n1 = sqrt(sum(a .* abs(xh).^2));
n2 = sqrt(sum(a .* abs(x).^2));
ip = max(abs(sum(a .* conj(xh) .* x)), abs(sum(a .* xh .* x)));
m.rho = ip / (n1*n2);
m.amp = n1 / n2;
m.relerr = sqrt(max(n1^2 + n2^2 - 2*ip, 0)) / n2;
end
