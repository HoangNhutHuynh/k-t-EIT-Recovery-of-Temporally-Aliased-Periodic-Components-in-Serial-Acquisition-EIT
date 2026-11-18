% run_E8_shape_baselines  E8: real data on thorax-shaped meshes with EIDORS temporal baselines.
% Inverse mesh warped to the shape_library outline (pig_23kg / neonate), anatomical lung ROI
% (lung_frac), regional nonlinearity index (nl_*: harmonic / fundamental amplitude, anterior vs
% posterior), and baselines ftc (reciprocity-completed), time_prior, kalman.
% Check fig_model.png of each run: electrode 1 must sit at the sternum, lungs inside the outline.
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;

rb = kt_config();
rb.rec.timing_source = 'estimate';
rb.rec.eidors_temporal = {'time_prior'};

cfg = rb; cfg.tag = 'E8_hfov_pigshape';
cfg.model.shape = 'pig_23kg';
cfg.real.files = {fullfile(pwd, 'Data', 'chb_hfov_porcine', 'CHB_pig06-55.get')};
cfg.real.cardiac_search = [1.5 3.5];
cfg.notes = 'HFOV on pig_23kg mesh + EIDORS baselines + NL index';
kt_exp_real(cfg);

cfg = rb; cfg.tag = 'E8_neonate_shape';
cfg.model.shape = 'neonate';
cfg.real.files = {fullfile(pwd, 'Data', 'if-neonate-ventilated', 'nm18.get')};
cfg.real.fs = 13; cfg.real.f0_search = [2.45 2.6]; cfg.real.n_harm = 3;
cfg.real.extra_excl = 0.4875*(1:12);
cfg.notes = 'neonate on neonate mesh + EIDORS baselines';
kt_exp_real(cfg);

files = {'p1130107', 'p1130122'};
for i = 1:numel(files)
    cfg = rb; cfg.tag = ['E8_' files{i} '_pigshape'];
    cfg.model.shape = 'pig_23kg';
    cfg.real.files = {fullfile(pwd, 'Data', 'pig_2003', [files{i} '.get'])};
    cfg.real.fs = 13; cfg.real.f0_search = [0.6 0.8]; cfg.real.n_harm = 12; cfg.real.block_sec = 10;
    cfg.notes = 'pig_2003 on pig_23kg mesh + EIDORS baselines';
    kt_exp_real(cfg);
end
