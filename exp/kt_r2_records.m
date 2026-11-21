function R = kt_r2_records(base)
%KT_R2_RECORDS The nine real records of the frame-removal experiment (Table IV).
%   R = kt_r2_records(base)   R(i).name, R(i).cfg (base with the record settings)
%   Settings as in run_FINAL_E7 (PEEP-trial piglet, HFOV pig, infant) and
%   run_E10_pig_torso (two piglets, Dargaville et al. 2010). g08-c035 is not
%   included: its cardiac line falls on the Nyquist frequency of the reduced rate.
%
% Author: Hoang Nhut Huynh

root = base.root;
R = struct('name', {}, 'cfg', {});

for f = {'p1130107', 'p1130122'}            % one newborn piglet, PEEP trial, after lavage / after surfactant
    c = base;
    c.real.files = {fullfile(root, 'Data', 'pig_2003', [f{1} '.get'])};
    c.real.fs = 13; c.real.f0_search = [0.6 0.8]; c.real.n_harm = 9;
    c.model.shape = 'pig_23kg';
    R(end+1) = struct('name', f{1}, 'cfg', c); %#ok<AGROW>
end

c = base;                                    % HFOV pig, 44 fps
c.real.files = {fullfile(root, 'Data', 'chb_hfov_porcine', 'CHB_pig06-55.get')};
c.real.fs = 44; c.real.f0_search = [8 12]; c.real.n_harm = 2; c.real.cardiac_search = [1.5 3.5];
c.model.shape = 'pig_23kg';
R(end+1) = struct('name', 'hfov', 'cfg', c);

c = base;                                    % ventilated infant, cardiac series
c.real.files = {fullfile(root, 'Data', 'if-neonate-ventilated', 'nm18.get')};
c.real.fs = 13; c.real.f0_search = [2.45 2.6]; c.real.n_harm = 2;
c.real.extra_excl = 0.4875*(1:12);
c.model.shape = 'neonate';
R(end+1) = struct('name', 'neonate', 'cfg', c);

% two piglets, PEEP trial 0 -> 25 -> 0 cmH2O in each record; cardiac line (or its 2nd harmonic) aliased
P = { 'g08-c003', [2.80 3.05], [5.75 5.95];
      'g08-c017', [2.45 2.70], [5.05 5.25];
      'g10-c005', [3.40 3.65], [];
      'g10-c019', [3.22 3.42], [];
      'g10-c029', [3.60 3.76], [] };
for r = 1:size(P, 1)
    c = base;
    c.real.files = {fullfile(root, 'Data', 'Pig_torso', [P{r, 1} '.get'])};
    c.real.fs = 13; c.real.f0_search = [0.30 0.42]; c.real.n_harm = 9;
    c.real.cardiac_search = P{r, 2};
    if ~isempty(P{r, 3}), c.real.extra_search = {P{r, 3}}; else, c.real.extra_search = {}; end
    c.model.shape = 'pig_23kg';
    R(end+1) = struct('name', P{r, 1}, 'cfg', c); %#ok<AGROW>
end
end
