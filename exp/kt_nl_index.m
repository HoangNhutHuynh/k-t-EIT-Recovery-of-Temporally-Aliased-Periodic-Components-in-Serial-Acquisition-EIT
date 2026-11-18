function rows = kt_nl_index(A, methods, names, C, I, roi, pre)
%KT_NL_INDEX Regional nonlinearity index: harmonic over fundamental oscillation amplitude.
%   rows = kt_nl_index(A, methods, names, C, I, roi, pre)
%   For each method, with |X| integrated over element area inside the ROI:
%     nl_total = sum_{h>=2} int|X_h| / int|X_1|
%     nl_hK    = int|X_K| / int|X_1|
%     nl_ant, nl_post : nl_total on the anterior / posterior half of the ROI,
%                       anterior = towards electrode 1 (sternum)
%     nl_ap    = nl_ant / nl_post
%   A harmonic placed at the wrong frequency or smeared by aliasing changes
%   this index, which is how a timing model can alter a physiological reading.
%   pre : '' (whole record) or 'blk' (per block).
%
% Author: Hoang Nhut Huynh

if ~isempty(pre), pre = [pre '_']; end
rows = cell(0, 6);
h1 = find(strcmp(names, 'h1'));
hs = find(~cellfun(@isempty, regexp(names, '^h\d+$')) & ~strcmp(names, 'h1') & [C.valid]);
if isempty(h1) || isempty(hs), return; end
p1 = mean(I.fmdl.nodes(I.fmdl.electrode(1).nodes, 1:2), 1);
u1 = p1 / norm(p1);
side = I.cent * u1.';
ant = roi & side > 0;
post = roi & side < 0;
for k = 1:numel(methods)
    m = methods{k};
    X = A.X.(m);
    if any(isnan(X(:, h1))), continue; end
    a = I.area(:);
    f1 = @(msk) sum(a(msk) .* abs(X(msk, h1)));
    fh = @(msk, h) sum(a(msk) .* abs(X(msk, h)));
    tot = 0; ta = 0; tp = 0;
    for h = hs
        if any(isnan(X(:, h))), continue; end
        r = fh(roi, h) / f1(roi);
        rows = [rows; {'', NaN, 'nl', m, [pre 'nl_' names{h}], r}]; %#ok<AGROW>
        tot = tot + r;
        ta = ta + fh(ant, h) / f1(ant);
        tp = tp + fh(post, h) / f1(post);
    end
    rows = [rows; {'', NaN, 'nl', m, [pre 'nl_total'], tot; '', NaN, 'nl', m, [pre 'nl_ant'], ta; ...
                   '', NaN, 'nl', m, [pre 'nl_post'], tp; '', NaN, 'nl', m, [pre 'nl_ap'], ta / tp}]; %#ok<AGROW>
end
end
