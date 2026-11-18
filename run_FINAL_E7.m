% run_FINAL_E7  Only the frame-rate reduction part of run_FINAL (key paper figure), about 10 min.
% Produces fig_decim_components.png: reference (full rate) vs naive, yerworth, ftc,
% time_prior and kt at the reduced rate, on thorax-shaped meshes.
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;
base = kt_config();
base.save_figs = true;
base.save_mat = true;
base.rec.eidors_temporal = {'time_prior'};

for f = {'p1130107', 'p1130122'}
    cfg = base; cfg.tag = ['FINAL_E7_' f{1}];
    cfg.real.files = {fullfile(pwd, 'Data', 'pig_2003', [f{1} '.get'])};
    cfg.real.fs = 13; cfg.real.f0_search = [0.6 0.8]; cfg.real.n_harm = 9;
    cfg.model.shape = 'pig_23kg';
    cfg.decim.fig_rows = {'h5', 'h6', 'h7', 'h8'};
    cfg.notes = 'final: pig_2003 13 -> 6.5 fps';
    kt_exp_decim(cfg);
end
cfg = base; cfg.tag = 'FINAL_E7_hfov';
cfg.real.files = {fullfile(pwd, 'Data', 'chb_hfov_porcine', 'CHB_pig06-55.get')};
cfg.real.n_harm = 2; cfg.real.cardiac_search = [1.5 3.5];
cfg.model.shape = 'pig_23kg';
cfg.notes = 'final: HFOV 44 -> 22 fps';
kt_exp_decim(cfg);
cfg = base; cfg.tag = 'FINAL_E7_neonate';
cfg.real.files = {fullfile(pwd, 'Data', 'if-neonate-ventilated', 'nm18.get')};
cfg.real.fs = 13; cfg.real.f0_search = [2.45 2.6]; cfg.real.n_harm = 2;
cfg.real.extra_excl = 0.4875*(1:12);
cfg.model.shape = 'neonate';
cfg.notes = 'final: neonate 13 -> 6.5 fps';
kt_exp_decim(cfg);
