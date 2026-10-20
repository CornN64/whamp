% Temporary equivalence probe (deleted after use).
root = 'c:\Users\WP\VSCODE\whamp';
addpath(fullfile(root, 'srcm'));

plasma.B = 100;
plasma.Species = { ...
    struct('m',16,'n',1,'t',10,'a',5,'b',0,'d',1,'vd',1), ...
    struct('m',0,'n',1,'t',100,'a',1,'b',0,'d',1,'vd',0)};

ip = struct('fstart',0.1,'kperp',[0 10 0],'kpar',[0.0022 10 0.0022], ...
    'varyKzFirst',0,'useLog',0,'maxIterations',50);

out = whamp_track_v2(plasma, ip);
fprintf('SIZE f = %s\n', mat2str(size(out.f)));
fprintf('f = %s\n', mat2str(out.f, 12));
fprintf('flagNoConvergence = %s\n', mat2str(out.flagNoConvergence));
fprintf('flagHeavyDamping  = %s\n', mat2str(out.flagHeavyDamping));
fprintf('certWinding = %s\n', mat2str(out.certWinding, 12));
fprintf('iterations  = %s\n', mat2str(out.iterations));
fprintf('kparUsed    = %s\n', mat2str(out.kparUsed, 12));
fprintf('kperpUsed   = %s\n', mat2str(out.kperpUsed, 12));
fprintf('DONE\n');
