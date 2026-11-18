function Mdl = kt_build_models(cfg, need_fwd)
%KT_BUILD_MODELS Coarse inverse model and (optionally) fine forward model.
%   Mdl = kt_build_models(cfg, need_fwd)
%   Mdl.inv / Mdl.fwd fields:
%     fmdl, vh (homogeneous voltages), J (normalised Jacobian, d(v/vh)/d sigma),
%     nodes, elems, cent, area, R (model radius), stim_idx, k_idx, n_stim
%   Mdl.RtR_laplace when cfg.rec.prior = 'laplace'.
%   cfg.model.shape ~= 'circle': inverse mesh warped to a shape_library outline
%   (kt_warp_shape), Mdl.inv.lung = anatomical lung mask (real data only; the
%   forward model used for simulation stays circular).
%   EIDORS calls: mk_common_model, mk_image, fwd_solve, calc_jacobian,
%   calc_RtR_prior, prior_laplace.
%
% Author: Hoang Nhut Huynh

n = cfg.acq.n_elec;
imdl = mk_common_model(cfg.model.inv, n);
lungs = {};
if ~strcmp(cfg.model.shape, 'circle') && ~need_fwd
    [imdl.fwd_model, lungs] = kt_warp_shape(imdl.fwd_model, cfg.model.shape, cfg.model.rotate_deg);
end
Mdl.imdl = imdl;
Mdl.inv = model_info(imdl.fwd_model, cfg.model.inv);
if ~isempty(lungs)        % anatomical lung mask on the inverse mesh
    Mdl.inv.lung = false(size(Mdl.inv.cent, 1), 1);
    for k = 1:numel(lungs)
        Mdl.inv.lung = Mdl.inv.lung | inpolygon(Mdl.inv.cent(:, 1), Mdl.inv.cent(:, 2), lungs{k}(:, 1), lungs{k}(:, 2));
    end
    Mdl.inv.lung_poly = lungs;
    fprintf('lung mask: %d of %d elements (%.0f%% of area)\n', nnz(Mdl.inv.lung), numel(Mdl.inv.lung), ...
        100*sum(Mdl.inv.area(Mdl.inv.lung))/sum(Mdl.inv.area));
end
Mdl.RtR_laplace = [];
if strcmp(cfg.rec.prior, 'laplace')
    imdl.RtR_prior = @prior_laplace;
    Mdl.RtR_laplace = calc_RtR_prior(imdl);
end
if need_fwd
    fm = mk_common_model(cfg.model.fwd, n);
    Mdl.fwd = model_info(fm.fwd_model, cfg.model.fwd);
    if ~isequal(Mdl.fwd.stim_idx, Mdl.inv.stim_idx) || ~isequal(Mdl.fwd.k_idx, Mdl.inv.k_idx)
        error('kt_build_models: forward and inverse models order the measurements differently');
    end
end
end

function I = model_info(fmdl, name)
fmdl.normalize_measurements = 0;   % normalisation is done here, not by EIDORS
img = mk_image(fmdl, 1);
vh = fwd_solve(img);
J = calc_jacobian(img);
I.fmdl = fmdl;
I.vh = vh.meas(:);
I.J = J ./ I.vh;
I.nodes = fmdl.nodes(:, 1:2);
I.elems = fmdl.elems;
if size(I.elems, 2) ~= 3
    error('kt_build_models: only 2-D triangular models are supported (got %d nodes/elem)', size(I.elems, 2));
end
a = I.nodes(I.elems(:, 1), :); b = I.nodes(I.elems(:, 2), :); c = I.nodes(I.elems(:, 3), :);
I.cent = (a + b + c) / 3;
I.area = abs((b(:,1)-a(:,1)).*(c(:,2)-a(:,2)) - (c(:,1)-a(:,1)).*(b(:,2)-a(:,2))) / 2;
I.R = max(sqrt(sum(I.nodes.^2, 2)));
[I.stim_idx, I.k_idx] = meas_map(fmdl);
I.n_stim = numel(fmdl.stimulation);
if numel(I.stim_idx) ~= numel(I.vh)
    error('kt_build_models: %d measurements in the stimulation struct but fwd_solve gave %d', ...
        numel(I.stim_idx), numel(I.vh));
end
fprintf('model %s: %d elements, %d measurements, %d injections\n', name, size(I.elems, 1), numel(I.vh), I.n_stim);
end

function [si, ki] = meas_map(fmdl)
% injection index and within-injection position of every measurement,
% in the order of fwd_solve output (EIDORS order: injection-major)
S = numel(fmdl.stimulation);
si = []; ki = [];
for s = 1:S
    r = size(fmdl.stimulation(s).meas_pattern, 1);
    si = [si; s*ones(r, 1)]; %#ok<AGROW>
    ki = [ki; (1:r).'];      %#ok<AGROW>
end
if isfield(fmdl, 'meas_select') && ~isempty(fmdl.meas_select) && numel(fmdl.meas_select) == numel(si)
    sel = logical(fmdl.meas_select(:));
    si = si(sel);
    ki = zeros(size(si));
    for s = 1:S
        ki(si == s) = 1:nnz(si == s);
    end
end
end
