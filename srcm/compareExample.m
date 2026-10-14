function compareExample()
%COMPAREEXAMPLE  Head-to-head over an adversarial line set.
%   Runs compareSolvers over nine sweep lines (three plasmas x three kperp)
%   and prints branch-hop and coverage counts for each solver.  Expected
%   result with the current code:
%       whamp         hops 13   missing 11
%       whamp_track   hops  0   missing 13
%   See compareSolvers.m for how a hop is detected.
addpath('c:\Users\WP\VSCODE\whamp\srcm');

% Three plasmas chosen to stress branch tracking: the golden drift case that
% exposes the bug, a strongly anisotropic one, and a fast-beam one.
golden = struct('B',100);
golden.Species = { ...
    struct('m',4,'n',1,'t',0.05,'a',5,'b',0,'d',1,'vd',1), ...
    struct('m',16,'n',0.5,'t',0.01,'a',1,'b',0,'d',1,'vd',0), ...
    struct('m',0,'n',1,'t',0.001,'a',1,'b',0,'d',1,'vd',0)};

aniso = struct('B',30);
aniso.Species = { ...
    struct('m',1,'n',0.6,'t',0.05,'a',2,'b',0,'d',1,'vd',0), ...
    struct('m',4,'n',0.4,'t',0.02,'a',3,'b',0,'d',1,'vd',0), ...
    struct('m',0,'n',1,'t',0.01,'a',1,'b',0,'d',1,'vd',0)};

beam = struct('B',100);
beam.Species = { ...
    struct('m',1,'n',0.7,'t',0.05,'a',1,'b',0,'d',1,'vd',2.5), ...
    struct('m',16,'n',0.3,'t',0.05,'a',1,'b',0,'d',1,'vd',0), ...
    struct('m',0,'n',1,'t',0.001,'a',1,'b',0,'d',1,'vd',0)};

plasmas = {golden, golden, golden, aniso, aniso, aniso, beam, beam, beam};
kperps  = {0.005, 0.01, 0.02, 0.005, 0.01, 0.02, 0.005, 0.01, 0.02};

lines = struct('plasma', {plasmas{:}}, 'kperp', {kperps{:}}, ...
    'z1', {0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08});

% NB: the extra cell level.  struct('solvers',{a,b}) would build a 1x2
% struct array, and opts.solvers would then expand to two arguments.
opts = struct('h', 0.008, 'refineFactor', 8, 'fstart', 0.4, ...
    'solvers', {{@whamp, @whamp_track, @whamp_global}}, 'verbose', false);
compareSolvers(lines, opts);
end
