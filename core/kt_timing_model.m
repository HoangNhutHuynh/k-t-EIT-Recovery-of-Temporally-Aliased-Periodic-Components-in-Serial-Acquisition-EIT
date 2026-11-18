function tau = kt_timing_model(stim_idx, k_idx, fs, P)
%KT_TIMING_MODEL Acquisition time offset tau_m (s) of every measurement in a frame.
%   tau = kt_timing_model(stim_idx, k_idx, fs, P)
%   stim_idx : M x 1 injection index (1..S) of each measurement
%   k_idx    : M x 1 position of the measurement inside its injection (1..K)
%   fs       : frame rate (Hz)
%   P.n_stim : S
%   P.mode   : 'per_stim' (all voltages of one injection read together)
%              'per_meas' (voltages read one after another)
%   P.duty   : fraction of the frame period spanned by the S injections (0..1)
%   P.shift  : injection slot that is read first (cyclic offset, 0..S-1)
%   P.dir    : +1 slots in increasing order, -1 reversed
%   Measurement m of frame f is taken at t = f/fs + tau(m), f = 0,1,...
%
% Author: Hoang Nhut Huynh

S = P.n_stim;
slot = mod(P.dir*(stim_idx(:) - 1) - P.shift, S);      % 0..S-1
switch P.mode
    case 'per_stim'
        u = slot / S;
    case 'per_meas'
        K = max(k_idx);
        u = (slot*K + (k_idx(:) - 1)) / (S*K);
    otherwise
        error('kt_timing_model: unknown mode "%s"', P.mode);
end
tau = P.duty * u / fs;
end
