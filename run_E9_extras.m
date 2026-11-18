% run_E9_extras  E9: robustness of the frame-rate reduction test (E7).
% E9a  regularisation: alpha 0.01 and 0.1 (E7 used 0.03), same alpha for every method
% E9b  stronger reduction: pig_2003 13 -> 4.33 fps (factor 3: h4-h9 aliased, ell = -1 and +1)
% E9c  h9 failure in E7 (folds to 0.18 Hz next to the 0.13 Hz slow line): add the slow line and
%      treat bins within 8 fs/F as one bin, so h9 and the line are solved jointly
%
% Author: Hoang Nhut Huynh

cd(fileparts(mfilename('fullpath')));
kt_addpaths;

files = {'p1130107', 'p1130122'};
hfov = fullfile(pwd, 'Data', 'chb_hfov_porcine', 'CHB_pig06-55.get');

% --- E9a: alpha
for a = [0.01 0.1]
    for i = 1:numel(files)
        cfg = pig_cfg(files{i});
        cfg.rec.alpha = a;
        cfg.tag = sprintf('E9a_%s_alpha%g', files{i}, a);
        cfg.notes = 'E7 with another alpha';
        kt_exp_decim(cfg);
    end
    cfg = kt_config();
    cfg.rec.eidors_temporal = {'time_prior'};
    cfg.real.files = {hfov};
    cfg.real.n_harm = 2; cfg.real.cardiac_search = [1.5 3.5];
    cfg.rec.alpha = a;
    cfg.tag = sprintf('E9a_hfov_alpha%g', a);
    cfg.notes = 'E7 HFOV with another alpha';
    kt_exp_decim(cfg);
end

% --- E9b: factor 3
for i = 1:numel(files)
    for tm = {'decimated', 'full'}
        cfg = pig_cfg(files{i});
        cfg.decim.factor = 3;
        cfg.decim.timing = tm{1};
        cfg.tag = sprintf('E9b_%s_D3_%s', files{i}, tm{1});
        cfg.notes = '13 -> 4.33 fps';
        kt_exp_decim(cfg);
    end
end

% --- E9c: slow line + near-collision joint solve
for i = 1:numel(files)
    cfg = pig_cfg(files{i});
    cfg.real.extra_search = {[0.08 0.2]};
    cfg.rec.tol_factor = 8;
    cfg.tag = sprintf('E9c_%s_slowline', files{i});
    cfg.notes = 'h9 + 0.13 Hz line solved jointly';
    kt_exp_decim(cfg);
end

function cfg = pig_cfg(f)
cfg = kt_config();
cfg.rec.eidors_temporal = {'time_prior'};
cfg.real.files = {fullfile(pwd, 'Data', 'pig_2003', [f '.get'])};
cfg.real.fs = 13; cfg.real.f0_search = [0.6 0.8]; cfg.real.n_harm = 9;
end
