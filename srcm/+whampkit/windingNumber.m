function n = windingNumber(plasma, p, z, x0, radius, numSamples)
%WINDINGNUMBER Number of zeros of D inside a circle, by the argument principle.
%   N = whampkit.windingNumber(PLASMA, P, Z, X0) counts the zeros of the
%   dispersion function inside the circle of radius RADIUS centred on X0,
%   counting poles negatively: N = Z - P.  It is the phase change of D
%   around the contour divided by 2*pi.
%
%   N = +1 therefore certifies a single simple zero, N = 0 certifies that
%   the interior contains neither a zero nor a pole, and N = -1 certifies a
%   pole.  Because the test is a property of D alone it is completely
%   scale-free, which matters here: |D| ranges over some thirty orders of
%   magnitude across the (p, z, x) domain, so no absolute residual
%   threshold can say whether a candidate point is a root.
%
%   N = whampkit.windingNumber(..., RADIUS) uses RADIUS instead of the
%   default 1e-4 * max(1, |X0|).
%
%   The contour is sampled ADAPTIVELY.  NUMSAMPLES (default 8) points are
%   used to start, and the count doubled until the largest phase step
%   between neighbouring samples falls below pi/2, so that no zero or pole
%   can hide in an unresolved arc.  A contour enclosing one simple zero has
%   an almost linear phase, so its first attempt already passes and costs
%   nine evaluations of D; a pole grazing the contour forces refinement up
%   to the NUMSAMPLES * 4 ceiling.  Capping the ceiling there means this is
%   never much dearer than the fixed 32-point contour it replaces, and
%   strictly more reliable, since a contour that coarse is now only accepted
%   once it has been shown to be resolved.  Passing NUMSAMPLES sets the
%   starting resolution; refinement from there still happens.
%
%   Points where the susceptibility expansion is refused (IERR) are treated
%   as "unknown" and the result is NaN, rather than silently counting as a
%   non-root.

if nargin < 5 || isempty(radius)
    radius = 1e-4 * max(1, abs(x0));
end
if nargin < 6 || isempty(numSamples)
    numSamples = 8;
end

% ORDER 1, not 2: the phase of D is all the argument principle uses, and
% dispersion() returns the same D either way while skipping the frequency
% derivative the Newton iteration needs and this test does not.
for refine = 0:2
    [n, worstStep] = local_contour(plasma, p, z, x0, radius, numSamples);
    if isnan(n) || worstStep <= pi / 2
        return
    end
    numSamples = 2 * numSamples;
end
end

function [n, worstStep] = local_contour(plasma, p, z, x0, radius, numSamples)
%CONTOUR  Winding of D around one closed polygon, and its coarsest arc.
%   Every sample shares one (p, z), so the whole polygon is evaluated in a
%   single vectorised pass through whampkit.evalDbatch rather than by
%   rebuilding the dispersion state once per corner.  The samples are the
%   vertices of a closed polygon, so the first is repeated at the end and
%   D is evaluated twice there; that keeps the phase increments between
%   neighbouring samples exactly what the scalar version produced.
theta = (0:numSamples) / numSamples * 2 * pi;
x = x0 + radius * complex(cos(theta), sin(theta));
[values, IERR] = whampkit.evalDbatch(plasma, p, z, x);
if any(IERR) || any(~isfinite(values)) || any(values == 0)
    n = NaN;
    worstStep = Inf;
    return
end

steps = angle(values(2:end) ./ values(1:end-1));
n = sum(steps) / (2 * pi);
n = round(n * 2) / 2;      % half-integers only if the contour is too coarse
worstStep = max(abs(steps));
end
