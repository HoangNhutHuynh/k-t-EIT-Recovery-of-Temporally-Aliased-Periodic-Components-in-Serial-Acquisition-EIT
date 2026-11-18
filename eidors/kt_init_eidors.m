function kt_init_eidors(cfg)
%KT_INIT_EIDORS Start EIDORS if it is not on the path yet.
%
% Author: Hoang Nhut Huynh

if exist('eidors_obj', 'file') ~= 2
    if isempty(cfg.eidors_startup) || ~exist(cfg.eidors_startup, 'file')
        error(['EIDORS is not on the path. Set cfg.eidors_startup to ' ...
               '<eidors folder>\startup.m in kt_config.m or in the run script.']);
    end
    run(cfg.eidors_startup);
end
fprintf('EIDORS: %s\n', which('eidors_obj'));
try
    eidors_msg('log_level', 1);
catch
end
end
