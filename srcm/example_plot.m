% Example: solve an anisotropic line with v2 and plot branch confidence.

aniso = struct('B', 30);
aniso.Species = { ...
    struct('m',1,'n',0.6,'t',0.05,'a',2,'b',0,'d',1,'vd',0), ...
    struct('m',4,'n',0.4,'t',0.02,'a',3,'b',0,'d',1,'vd',0), ...
    struct('m',0,'n',1, 't',0.01,'a',1,'b',0,'d',1,'vd',0)};

ip = struct('fstart',0.4,'kperp',0.005,'kpar',[0 0.008 0.08], ...
            'maxIterations',50);

out = whamp_track(aniso, ip);

fprintf('found %d   suspect %d   notCertified %d   certUnknown %d\n', ...
    nnz(out.flagSolutionFound), nnz(out.flagSuspectBranch), ...
    nnz(out.flagNotCertified), nnz(out.flagCertUnknown));

figure;
P = whamp_track_plot(out, 'ShowPredicted', true);
title('aniso k_{perp}=0.005  —  red X = suspect branch');

saveas(gcf, 'example_plot.png');
fprintf('handles: %s\n', strjoin(fieldnames(P), ', '));
