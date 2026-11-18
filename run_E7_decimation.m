% run_E7_decimation  E7: frame-rate reduction on real data (ground truth from the full rate).
% Harmonics below Nyquist at the recorded rate become aliased after keeping every 2nd frame;
% their full-rate images are the reference. Two timing variants per record:
% schedule re-estimated at the low rate ('decimated', realistic) and taken from the full rate ('full').
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;

base = kt_config();
base.rec.eidors_temporal = {'time_prior'};

% pig_2003, 13 -> 6.5 fps: h5-h9 (3.5-6.3 Hz) become aliased
for f = {'p1130107', 'p1130122'}
    for tm = {'decimated', 'full'}
        cfg = base;
        cfg.real.files = {fullfile(pwd, 'Data', 'pig_2003', [f{1} '.get'])};
        cfg.real.fs = 13; cfg.real.f0_search = [0.6 0.8]; cfg.real.n_harm = 9;
        cfg.decim.timing = tm{1};
        cfg.tag = sprintf('E7_%s_%s', f{1}, tm{1});
        cfg.notes = 'pig_2003 13 -> 6.5 fps';
        kt_exp_decim(cfg);
    end
end

% HFOV pig, 44 -> 22 fps: h2 (20.6 Hz) becomes aliased; cardiac line kept
for tm = {'decimated', 'full'}
    cfg = base;
    cfg.real.files = {fullfile(pwd, 'Data', 'chb_hfov_porcine', 'CHB_pig06-55.get')};
    cfg.real.n_harm = 2; cfg.real.cardiac_search = [1.5 3.5];
    cfg.decim.timing = tm{1};
    cfg.tag = sprintf('E7_hfov_%s', tm{1});
    cfg.notes = 'HFOV 44 -> 22 fps';
    kt_exp_decim(cfg);
end

% neonate, 13 -> 6.5 fps: cardiac h2 (4.9 Hz) becomes aliased
for tm = {'decimated', 'full'}
    cfg = base;
    cfg.real.files = {fullfile(pwd, 'Data', 'if-neonate-ventilated', 'nm18.get')};
    cfg.real.fs = 13; cfg.real.f0_search = [2.45 2.6]; cfg.real.n_harm = 2;
    cfg.real.extra_excl = 0.4875*(1:12);
    cfg.decim.timing = tm{1};
    cfg.tag = sprintf('E7_neonate_%s', tm{1});
    cfg.notes = 'neonate 13 -> 6.5 fps';
    kt_exp_decim(cfg);
end
