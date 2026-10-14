function [roots, info] = findAllRoots(plasma, p, z, box, options)
%FINDALLROOTS Enumerate the dispersion roots in a frequency window.
%   ROOTS = whampkit.findAllRoots(PLASMA, P, Z, BOX) returns the distinct
%   roots of D(z, w) inside BOX = [wmin wmax imin imax], ordered with ESORT.
%
%   [ROOTS, INFO] also reports:
%     .count       number of distinct roots
%     .evals       dispersion evaluations used
%     .seeds       seed samples tried
%     .stable      true when doubling the seed density found nothing new
%     .unstable    roots that appeared only at the finer density
%
%   Why this exists: at a fixed (kperp, kpar) the WHAMP dispersion function
%   generally has several roots at once, and a single-seed Newton solver
%   reports whichever one its starting value happens to lie in the basin of.
%   Enumerating the whole set first turns "which root is the continuation of
%   the previous point?" into a well-posed matching problem, which is what
%   whamp_global uses to stop the reported branch from jumping.
%
%   Method
%     Multi-start Newton from a dense set of seeds laid along a line inside
%     the box, followed by deduplication.  Two details make it work, and
%     both are consequences of this particular D:
%
%     * |D| cannot be thresholded.  The dispersion function is written in a
%       form whose magnitude runs 1e20..1e40 across the band of interest --
%       even the reference solver's own converged root sits near |D| ~ 1e35.
%       So Newton's stopping test is on the *correction*, and acceptance is
%       decided by the residue test below rather than by |D|.
%
%     * The roots lie essentially on the real axis, interleaved with the
%       real-axis gyroharmonic poles.  A scan along a line slightly off the
%       axis resolves them; a scan on the axis cannot, because the |D| dip
%       at a root is ~1e-12 wide.
%
%     * Newton is step-capped.  D is meromorphic with a pole at every
%       species gyroharmonic (bl/y with y = x - dl*z, and cot(pi*y) inside
%       r_asymptotic).  An uncapped step started near a pole flies across
%       the plane and lands in an arbitrary basin, which makes the
%       enumeration depend on seed order and, at high seed density,
%       manufactures spurious far-away "roots".
%
%   Completeness.  A contour (argument-principle) count is not used: the
%   gyroharmonic poles sit on the real axis inside the region of interest,
%   so no rectangle can be drawn there without either crossing a pole or
%   excluding the roots that live between it.  Instead the scan is repeated
%   at double the seed density and the two result sets compared -- a root
%   that survives refinement is stable, and INFO.stable says whether the
%   coarser pass was already complete.
%
%   Options (struct, all optional):
%     .seeds      Number of seeds on the coarse pass.     Default 120
%     .tol        Relative Newton tolerance.              Default 1e-10
%     .maxIter    Newton iteration cap.                   Default 60
%     .dedupTol   Two roots closer than this are one.     Default 1e-6
%     .stepCap    Max Newton step, fraction of |w|.       Default 0.1
%     .refine     Do the density-doubling check.          Default true
%     .extraSeeds Additional seed frequencies (row vector), always tried.

if nargin < 5 || isempty(options)
    options = struct();
end
nSeeds = local_get(options, 'seeds', 120);
tol = local_get(options, 'tol', 1e-10);
maxIter = local_get(options, 'maxIter', 60);
dedupTol = local_get(options, 'dedupTol', 1e-6);
stepCap = local_get(options, 'stepCap', 0.1);
doRefine = local_get(options, 'refine', true);
extraSeeds = local_get(options, 'extraSeeds', []);

box = local_sanitiseBox(plasma, box, z);
wmin = box(1);
wmax = box(2);

% Scan line.  The roots sit essentially on the real axis (imaginary parts
% 1e-14 .. 1e-4) interleaved with real-axis gyroharmonic poles, and the
% |D| dip at a root is far too narrow to resolve by sampling the axis
% itself.  What works is seeding ON the axis and letting Newton's complex
% corrections walk off it to the root: the seed only has to be in the right
% basin.  Seeds on the axis are also always inside the susceptibility's
% validity cone.
imLine = local_get(options, 'scanImag', 0);

[coarse, evals, samples] = local_scan(plasma, p, z, wmin, wmax, imLine, ...
    nSeeds, tol, maxIter, stepCap, dedupTol, extraSeeds);

stable = true;
unstable = zeros(1, 0);
if doRefine
    [fine, evals2, samples2] = local_scan(plasma, p, z, wmin, wmax, imLine, ...
        2 * nSeeds, tol, maxIter, stepCap, dedupTol, extraSeeds);
    evals = evals + evals2;
    samples = samples + samples2;
    for k = 1:numel(fine)
        if ~any(abs(coarse - fine(k)) <= dedupTol * max(1, abs(fine(k))))
            unstable(end+1) = fine(k); %#ok<AGROW>
        end
    end
    stable = isempty(unstable);
    roots = local_esort(fine);
else
    roots = local_esort(coarse);
end

info = struct('count',numel(roots),'evals',evals,'seeds',samples, ...
    'stable',stable,'unstable',unstable);
end

% ------------------------------------------------------------------ helpers

function [roots, evals, samples] = local_scan(plasma, p, z, wmin, wmax, ...
        imLine, n, tol, maxIter, stepCap, dedupTol, extraSeeds)
%SCAN Multi-start Newton along the horizontal scan line.
seeds = linspace(wmin, wmax, n) + 1i * imLine;
if ~isempty(extraSeeds)
    seeds = [seeds, extraSeeds(:).' + 1i * imLine];
end

roots = zeros(1, 0);
evals = 0;
samples = numel(seeds);
for k = 1:samples
    [root, converged, cost] = local_newton(plasma, p, z, seeds(k), ...
        tol, maxIter, stepCap);
    evals = evals + cost;
    if ~converged
        continue
    end
    if real(root) < wmin || real(root) > wmax
        continue
    end
    if any(abs(roots - root) <= dedupTol * max(1, abs(root)))
        continue
    end
    roots(end+1) = root; %#ok<AGROW>
end
end

function [root, converged, evals] = local_newton(plasma, p, z, x0, tol, maxIter, stepCap)
%NEWTON Step-capped Newton on D(w), accepted by the residue test.
%   The step cap is the important part.  Newton applied to a meromorphic
%   function is attracted to zeros, but a seed sitting close to a pole
%   produces a correction that leaps over the whole plane; capping keeps the
%   iterate in the neighbourhood of the singularity structure it started in,
%   so the set of roots found depends on where the seeds are and not on
%   numerical luck.
root = x0;
converged = false;
evals = 0;
for it = 1:maxIter
    d = whampkit.evalD(plasma, p, z, root, 2);
    evals = evals + 1;
    if d.IERR || ~isfinite(d.D) || ~isfinite(d.DX) || d.DX == 0
        return
    end
    correction = d.D / d.DX;
    if ~isfinite(correction)
        return
    end
    cap = stepCap * max(abs(root), 1e-6);
    if abs(correction) > cap
        correction = correction / abs(correction) * cap;
    end
    root = root - correction;
    if abs(correction) <= tol * max(abs(root), 1e-12)
        break
    end
end
converged = local_isZero(plasma, p, z, root);
evals = evals + 2;
end

function ok = local_isZero(plasma, p, z, root)
%ISZERO  Scale-free test that ROOT is a zero of D rather than a pole.
%   For D ~ a(w - x0)^m the logarithmic derivative satisfies
%   (w - x0) D'(w) / D(w) -> m, a positive integer at a zero of order m and
%   -1 at a simple pole.  Probing at a radius far larger than the Newton
%   residual but far smaller than the spacing to the next singularity gives
%   an accurate, unit-free verdict -- which is the only kind available, since
%   |D| itself is of order 1e35 even at a genuine root.
probe = 1e-6 * max(1, abs(root));
d = whampkit.evalD(plasma, p, z, root + probe, 2);
if d.IERR || ~isfinite(d.D) || ~isfinite(d.DX) || d.D == 0
    ok = false;
    return
end
residue = probe * d.DX / d.D;
ok = abs(residue - round(real(residue))) < 0.05 && round(real(residue)) >= 1;
end

function box = local_sanitiseBox(plasma, box, z)
%SANITISEBOX Clip the search box to where the susceptibility is defined.
%   whampkit.susceptibility raises IERR once imag(x) drops below
%   -0.6*|real(x) - z*vd|, so a box reaching far below the real axis spends
%   most of its seeds on evaluations that are rejected outright.  Growth
%   rates in this model are tiny next to the real frequency, so clipping
%   costs nothing in practice.
imin = box(3);
if imin >= 0
    return
end
vd = 0;
if isfield(plasma, 'drift') && ~isempty(plasma.drift)
    vd = max(abs(plasma.drift(:)));
end
rx = min(abs([box(1) box(2)] - z * vd));
if box(1) <= z * vd && z * vd <= box(2)
    rx = 0;
end
limit = -0.6 * rx;
if imin < limit
    box(3) = 0.9 * limit;
end
end

function r = local_esort(x)
%ESORT Order a root set the way the Control System Toolbox does.
%   esort sorts by imaginary part and then real part, giving a listing that
%   is independent of the order the roots happened to be discovered in.  That
%   reproducibility is what lets two neighbouring k points be compared
%   index by index.  Falls back to an equivalent local sort if unavailable.
if isempty(x)
    r = x;
    return
end
if exist('esort', 'file') == 2 || exist('esort', 'builtin') == 2
    try
        r = esort(x);
        return
    catch
    end
end
[im, order] = sort(imag(x(:).'));
re = real(x(order));
[~, order2] = sort(2^52 * im + re);
r = x(order(order2));
end

function v = local_get(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end
