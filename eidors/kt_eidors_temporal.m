function X = kt_eidors_temporal(imdl, Iinv, mrow, dv, fs, bins, win, method, alpha)
%KT_EIDORS_TEMPORAL Baselines from EIDORS temporal solvers, projected on the component bins.
%   X = kt_eidors_temporal(imdl, Iinv, mrow, dv, fs, bins, win, method, alpha)
%   method 'time_prior' : inv_solve_time_prior (time_steps 2, weight 0.8)
%          'kalman'     : inv_solve_diff_kalman
%   Frame images (no within-frame timing model) are reconstructed from the
%   stored rows mrow of the full model (Iinv = Mdl.inv before restriction) and
%   their windowed DTFT at bins gives N x numel(bins) component images.
%   Tikhonov prior, hyperparameter from the same alpha rule as kt_recon_matrix
%   (unweighted). Scale differs from the k-t images: compare with rho, not amp.
%   Returns [] and prints the reason on any EIDORS error.
%
% Author: Hoang Nhut Huynh

X = [];
try
    fm = imdl.fwd_model;
    nsel = numel(Iinv.vh);
    if isfield(fm, 'meas_select') && ~isempty(fm.meas_select) && numel(fm.meas_select) ~= nsel
        full = find(fm.meas_select);
        mask = false(size(fm.meas_select));
        mask(full(mrow)) = true;
    else
        mask = false(nsel, 1);
        mask(mrow) = true;
    end
    [srt, ord] = sort(mrow(:));
    fm.meas_select = mask;
    fm.normalize_measurements = 1;
    imdl.fwd_model = fm;
    imdl.reconst_type = 'difference';
    imdl.jacobian_bkgnd.value = 1;
    imdl.RtR_prior = @prior_tikhonov;
    Jn = Iinv.J(srt, :);
    imdl.hyperparameter.value = sqrt(alpha * trace(Jn.' * Jn) / size(Jn, 2));
    switch method
        case 'time_prior'
            imdl.solve = @inv_solve_time_prior;
            imdl.inv_solve_time_prior.time_steps = 2;
            imdl.inv_solve_time_prior.time_weight = 0.8;
        case 'kalman'
            imdl.solve = @inv_solve_diff_kalman;
        otherwise
            error('unknown method %s', method);
    end
    vh = Iinv.vh(srt);
    vi = vh .* (1 + dv(ord, :));
    img = inv_solve(imdl, vh, vi);
    X = kt_spectrum(img.elem_data, fs, bins, win);
    fprintf('EIDORS %s baseline: %d frames reconstructed\n', method, size(img.elem_data, 2));
catch err
    fprintf('EIDORS %s baseline skipped: %s\n', method, err.message);
end
end
