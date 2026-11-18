% run_E1_sim  E1: main simulation (HFOV 10 Hz + 20/30 Hz harmonics + cardiac, fs = 44 Hz).
% Each block below is one logged run. Copy a block, change tag/notes/fields to try a variant.
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;

% E1a: oracle timing (tests the data model itself)
cfg = kt_config();
cfg.tag = 'E1a_oracle';
cfg.notes = 'true acquisition times; kt should equal simref';
kt_exp_sim(cfg);

% E1b: timing estimated from the coherence at f0 (what real data allows)
cfg = kt_config();
cfg.tag = 'E1b_est';
cfg.rec.timing_source = 'estimate';
cfg.notes = 'timing from 10 Hz coherence';
kt_exp_sim(cfg);

% E1c: nonlinear forward model (one fwd_solve per injection; slower, a few minutes)
cfg = kt_config();
cfg.tag = 'E1c_nonlin';
cfg.sim.forward = 'nonlinear';
cfg.rec.timing_source = 'estimate';
cfg.notes = 'nonlinear forward, estimated timing';
kt_exp_sim(cfg);
