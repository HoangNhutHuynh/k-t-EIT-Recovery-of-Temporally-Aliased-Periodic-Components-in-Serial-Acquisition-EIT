function [w, s2, probe] = kt_noise_floor(dv, fs, f_excl, win, guard)
%KT_NOISE_FLOOR Per-channel noise variance from spectral bins without signal.
%   [w, s2] = kt_noise_floor(dv, fs, f_excl, win, guard)
%   Probes frequencies between 0.5 Hz and fs/2 - 0.5 Hz that are more than
%   guard (Hz) away from every frequency in f_excl, and takes the median
%   periodogram value per channel. For complex Gaussian noise the mean of
%   |Y|^2 equals median/ln(2). w = 1./s2 are the data weights (W = diag(w)),
%   so the data term of kt_solve is in chi-square units.
%
% Author: Hoang Nhut Huynh

F = size(dv, 2);
df = fs/F;
probe = 0.5:df:(fs/2 - 0.5);
keep = true(size(probe));
for k = 1:numel(f_excl)
    keep = keep & abs(probe - f_excl(k)) > guard;
end
probe = probe(keep);
if numel(probe) < 10
    error('kt_noise_floor: only %d probe bins left (record too short or guard too wide)', numel(probe));
end
Y = kt_spectrum(dv, fs, probe, win);
s2 = median(abs(Y).^2, 2) / log(2);
w = 1 ./ s2;
end
