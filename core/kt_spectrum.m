function Y = kt_spectrum(dv, fs, f, win)
%KT_SPECTRUM Windowed DTFT of every measurement channel at chosen frequencies.
%   Y = kt_spectrum(dv, fs, f, win)
%   dv  : M x F data (channel mean is removed here)
%   f   : frequencies (Hz), need not lie on FFT bins
%   win : 'hann' (default) or 'rect'
%   Y   : M x numel(f). Scaled so that a channel a*cos(2*pi*f*t + p),
%         t = frame start times, gives Y = a*exp(1i*p).
%
% Author: Hoang Nhut Huynh

if nargin < 4, win = 'hann'; end
F = size(dv, 2);
n = (0:F-1).';
switch win
    case 'hann', w = 0.5 - 0.5*cos(2*pi*(n + 0.5)/F);
    case 'rect', w = ones(F, 1);
    otherwise, error('kt_spectrum: unknown window "%s"', win);
end
dv = dv - mean(dv, 2);
E = exp(-1i*2*pi*(n/fs)*f(:).');          % F x K
Y = 2 * dv * (w .* E) / sum(w);
end
