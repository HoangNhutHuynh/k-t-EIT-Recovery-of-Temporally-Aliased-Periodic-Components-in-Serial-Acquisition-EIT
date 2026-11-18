% run_E6_pig2003  E6: third real record set, pig_2003 (Goe-MF II, 13 frames/s, 208 measurements, 115 s).
% A periodic component at f0 = 0.70 Hz with a long harmonic series (cardiac or ventilator; checked on
% the spectrum: harmonics 1-9 unaliased, 10th (7.0 Hz -> 6.0 Hz, ell = -1) about 6-7x above the local
% floor, 12th (8.4 Hz -> 4.6 Hz) 2-4x, 11th and 13th marginal). Weak aliased lines: a hard test.
% 10 s blocks: with 5 s blocks the bins of h10 (5.99 Hz) and h8 (5.61 Hz) would be merged.
% Runtime: free timing + CV with 12 components, a few minutes per file.
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;

files = {'p1130107.get', 'p1130122.get'};
for i = 1:numel(files)
    cfg = kt_config();
    cfg.real.files = {fullfile(pwd, 'Data', 'pig_2003', files{i})};
    cfg.real.fs = 13;
    cfg.real.f0_search = [0.6 0.8];
    cfg.real.n_harm = 12;
    cfg.real.block_sec = 10;
    cfg.rec.timing_source = 'estimate';
    cfg.rec.timing_method = 'free';
    cfg.real.cv_timing = true;
    cfg.tag = ['E6_pig2003_' strtok(files{i}, '.')];
    cfg.notes = 'third subject: f0 0.70 Hz series, h10-h12 aliased (weak)';
    kt_exp_real(cfg);
end
