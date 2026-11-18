% run_E1d_phase  E1d: lungs oscillating out of phase (pendelluft-like time lag).
% The timing estimator assumes an in-phase oscillation at f0; this shows how
% coherence, timing error and the harmonic images degrade when that is false.
% A time lag dt gives phase h*2*pi*f0*dt at harmonic h, so lag at f0 = L deg
% means 2L at 20 Hz and 3L at 30 Hz.
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;

for L = [15 30 60 90]
    cfg = kt_config();
    cfg.tag = sprintf('E1d_lag%d', L);
    cfg.rec.timing_source = 'estimate';
    for c = 1:numel(cfg.sim.components)
        h = cfg.sim.components(c).harm;
        if h > 0, cfg.sim.components(c).phase_lr_deg = h*L; end
    end
    cfg.notes = sprintf('left-right lag %d deg at f0, estimated timing', L);
    kt_exp_sim(cfg);
end
