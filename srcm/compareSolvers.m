function compareSolvers(lines, opts)
%COMPARESOLVERS  Measure branch fidelity of the three WHAMP solvers.
%
%   How a solver is judged
%   ----------------------
%   Ground truth for "which root continues this branch" cannot come from a
%   global root finder here: in cold multi-species plasmas the susceptibility
%   expansion fails over most of the complex plane, so an enumerator cannot
%   sample the box reliably.
%
%   It does not need to.  The baseline loses a branch only because its
%   one-sided linear extrapolation seed -- built from the previous point --
%   falls into another basin.  Shrink the step and that seed becomes
%   arbitrarily accurate, so the baseline on a step N times finer tracks the
%   branch the coarse run should have followed.  It costs milliseconds, needs
%   no extra assumptions, and is the same code the user runs.
%
%   The comparison is then essentially exact.  Both runs report converged
%   roots of the same equation at the *same* kpar values (the fine grid
%   contains the coarse grid), so when they agree on the branch they agree to
%   ~1e-12, and when the coarse run hops the discrepancy is of order the
%   spacing between branches.  There is nothing to tune in between.
%
%   compareSolvers(LINES)     LINES: struct array with .plasma .kperp .z1
%   compareSolvers(LINES,OPT) OPT.h (coarse step, default 0.004),
%                             OPT.refineFactor (default 8), OPT.fstart,
%                             OPT.solvers, OPT.verbose, OPT.tol
%
%   Also see: whamp_track.m, whamp_global.m

addpath('c:\Users\WP\VSCODE\whamp\srcm');

if nargin < 2, opts = struct(); end
N = local_get(opts, 'refineFactor', 8);
h = local_get(opts, 'h', 0.004);
fstart = local_get(opts, 'fstart', 0.1);
solvers = local_get(opts, 'solvers', {@whamp, @whamp_track, @whamp_global});
verbose = local_get(opts, 'verbose', true);
tol = local_get(opts, 'tol', 1e-6);

tot = struct('name',{{}}, 'hop',{}, 'miss',{}, 'checked',{}, 'skipped',{});
for s = 1:numel(solvers)
    tot(s).name = func2str(solvers{s});
    tot(s).hop = 0;
    tot(s).miss = 0;
    tot(s).checked = 0;
    tot(s).skipped = 0;
end

for li = 1:numel(lines)
    L = lines(li);
    fprintf('\n%s\nLINE %d  kperp=%g  kpar 0:%g:%g\n%s\n', ...
        repmat('=', 1, 76), li, L.kperp, h, L.z1, repmat('=', 1, 76));

    fine = run_solver(@whamp, L, fstart, h/N, L.z1);
    if numel(fine.fv) < 5
        fprintf('  reference too short (%d points), skipping\n', numel(fine.fv));
        continue
    end
    % The reference must itself be a smooth single branch, or it proves nothing.
    fprintf('  reference: %d/%d converged, max extrapolation residual %.1e\n', ...
        numel(fine.fv), fine.nPoints, local_roughness(fine.fv));

    for s = 1:numel(solvers)
        R = run_solver(solvers{s}, L, fstart, h, L.z1);
        [hop, miss, worst, checked, skipped, msg] = judge(R, fine, tol);
        tot(s).hop = tot(s).hop + hop;
        tot(s).miss = tot(s).miss + miss;
        tot(s).checked = tot(s).checked + checked;
        tot(s).skipped = tot(s).skipped + skipped;
        % 'checked' is how many points the reference could actually judge.  A
        % solver can only be scored on those; if the reference died early, a
        % hop further along the line is invisible here and the hop count is a
        % lower bound.  That is why 'checked' is printed alongside.
        fprintf('  %-13s hops %d  missing %d  worst dev %.2e  %5.1fs  (%d checked, %d unjudgeable)\n', ...
            tot(s).name, hop, miss, worst, R.elapsedSecs, checked, skipped);
        if verbose
            for k = 1:numel(msg)
                fprintf('      %s\n', msg{k});
            end
        end
    end
end

fprintf('\n%s\nTOTALS over %d lines\n%s\n', repmat('=', 1, 76), numel(lines), ...
    repmat('=', 1, 76));
for s = 1:numel(tot)
    fprintf('  %-13s hops %d   missing %d\n', tot(s).name, tot(s).hop, tot(s).miss);
end
end

function v = local_get(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end

function r = local_roughness(fv)
%LOCALROUGHNESS  Max |f(i) - 2f(i-1) + f(i-2)|: ~0 for a smooth branch.
if numel(fv) < 3
    r = 0;
else
    r = max(abs(fv(3:end) - 2*fv(2:end-1) + fv(1:end-2)));
end
end

function R = run_solver(fh, L, fstart, h, z1)
P = struct('fstart', fstart, 'kperp', L.kperp, 'kpar', [0 h z1], ...
    'varyKzFirst', 0, 'useLog', 0, 'maxIterations', 50);
t = tic;
R = fh(L.plasma, P);
R.elapsedSecs = toc(t);
ok = R.flagSolutionFound(:) ~= 0;
f = R.f(:);
k = R.kpar(:);
R.nPoints = numel(f);
keep = ok & isfinite(f) & (abs(f) > 0);
R.k = k(keep);
R.fv = f(keep);
end

function [hop, miss, worst, checked, skipped, msg] = judge(R, fine, tol)
%JUDGE  Points where the run is not on the reference branch.
%
%   Matched at identical kpar (the refined grid contains the coarse one), so
%   both sides are converged roots of the same dispersion relation and any
%   disagreement beyond round-off means a different root.
%
%   A point the reference did not reach cannot be judged and is counted in
%   SKIPPED rather than silently ignored: without that, a solver that hops
%   beyond where the reference gave out would score zero hops, which is how
%   this harness first understated the baseline.
msg = {};
hop = 0;
worst = 0;
checked = 0;
skipped = 0;
miss = R.nPoints - numel(R.k);

for i = 1:numel(R.fv)
    [~, j] = min(abs(fine.k - R.k(i)));
    if ~(abs(fine.k(j) - R.k(i)) <= 1e-9 * max(1, abs(R.k(i))))
        skipped = skipped + 1;
        continue
    end
    checked = checked + 1;
    dev = abs(R.fv(i) - fine.fv(j));
    worst = max(worst, dev);
    if dev > tol * max(1, abs(fine.fv(j)))
        hop = hop + 1;
        msg{end+1} = sprintf('kpar=%-9g f=%-14.7g  branch=%-14.7g  dev %.2e', ...
            R.k(i), real(R.fv(i)), real(fine.fv(j)), dev); %#ok<AGROW>
    end
end
end
